#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Module 5's long-flagged remaining item: GMRES/BiCGSTAB. Its concrete
// unlock trigger, identified during a planning checkpoint but not acted on
// until now: every linear system solved anywhere in this project so far
// (pressure projection, Poiseuille/heat diffusion, k/omega transport's own
// Patankar-linearized sweeps) is symmetric, or solved by plain Gauss-Seidel
// which doesn't care -- so Conjugate Gradient has always been the right,
// sufficient tool, and no genuinely non-symmetric matrix has ever needed to
// be assembled. Steady 1D convection-diffusion with a first-order upwind
// convective term is the simplest case that actually needs one: upwinding
// weights the two neighbor coefficients unevenly depending on flow
// direction, breaking the symmetry every purely-diffusive stencil here has.
//
// Governing equation: u * dphi/dx = nu * d2phi/dx2 + source, steady,
// Dirichlet at both ends (phi(0)=leftValue, phi(L)=rightValue). Finite
// volume discretization: central differencing for the diffusive face flux
// (same as every diffusion solver here), first-order upwind for the
// convective face flux -- the standard choice that stays bounded/stable at
// any cell Peclet number, unlike central differencing for convection
// (which needs Pe_cell < 2), at the cost of numerical ("false") diffusion
// measured in the tests. Same ghost-mirror Dirichlet wall convention as
// MixingLengthChannelFlowSolver1D (the wall sits exactly half a cell beyond
// the boundary cell center).
//
// **Exact solution, derived here (not recalled) from the governing ODE**,
// for the constant-coefficient, source-free case:
// u*phi' = nu*phi'' has general solution phi = A*exp(u*x/nu) + B; applying
// the two Dirichlet conditions gives
//   phi(x) = leftValue + (rightValue-leftValue) * (exp(u*x/nu)-1) / (exp(u*L/nu)-1)
// (an ordinary linear-ODE boundary-value solution, not a benchmark table --
// consistent with this project's practice of deriving validation targets
// rather than reciting literature). First-order upwind's own false
// diffusion means the discrete solution only matches this closely at
// low-to-moderate cell Peclet number; how closely was measured directly,
// not guessed, before picking a test tolerance.
//
// **Spatially varying coefficients** (setVelocityField/setDiffusivityField)
// are supported, and are not a contrived generalization: every turbulence
// closure in this project computes an effective diffusivity nu + nu_t that
// varies cell to cell, so a variable-coefficient convection-diffusion
// operator is exactly the shape a real implicit turbulent-transport solve
// would take. They also matter for a specific numerical reason -- see
// setPreconditioner() below.
class ImplicitConvectionDiffusionSolver1D {
public:
    // Preconditioner applied (on the left) to the Krylov solvers below:
    // they solve M^-1 A x = M^-1 b instead of A x = b, which converges
    // faster when M approximates A but is much cheaper to invert.
    //
    // **Jacobi is deliberately useless on a *constant*-coefficient problem,
    // and that is worth stating rather than discovering later.** This
    // system's diagonal coefficient is |u| + 2*nu/h at every interior cell;
    // with constant u and nu that is the same number in every row (measured:
    // exactly 42.0 for the default test case). M = D is then a scalar
    // multiple of the identity, so M^-1 A is just A rescaled -- and Krylov
    // methods are invariant under scalar rescaling (identical Krylov
    // subspace, identical minimizer). Jacobi only becomes a real
    // preconditioner once the diagonal actually varies, which is what
    // setVelocityField()/setDiffusivityField() make possible.
    //
    // **IncompleteLU (ILU(0)) is exact for this operator, and that is a
    // property of 1D, not a general claim.** ILU(0) keeps only the fill-in
    // pattern of A itself; a 1D finite-volume stencil is tridiagonal, and
    // LU factorization of a tridiagonal matrix produces no fill-in outside
    // that pattern, so nothing is dropped and ILU(0) *is* the exact LU
    // factorization. M^-1 A is then exactly the identity and a Krylov
    // solver converges in a single iteration. In 2D/3D the stencil is no
    // longer tridiagonal, ILU(0) does drop real fill-in, and it becomes a
    // genuine approximation -- so the one-iteration result measured here
    // must not be read as "ILU(0) always solves it instantly".
    enum class Preconditioner { None, Jacobi, IncompleteLU };

    // How the convected value is reconstructed at a face.
    //
    // **The assembled matrix is upwind whichever of these is chosen**, and
    // the difference is carried on the right-hand side by deferred
    // correction, re-evaluated over outer sweeps. That is the same mechanism
    // UnstructuredFvmBase uses for its non-orthogonal term, for the same two
    // reasons: the matrix keeps the M-matrix property that makes the Krylov
    // solve well behaved, and a *nonlinear* scheme (any limiter) cannot be
    // put in a linear operator at all.
    enum class ConvectionScheme {
        // phi_f = phi_upwind. Bounded at any cell Peclet number, first-order
        // accurate, and its error is false diffusion -- which is why this
        // class's own tests measure that error rather than assuming it.
        FirstOrderUpwind,
        // phi_f = (phi_C + phi_D)/2. Second-order accurate and
        // **unconditionally oscillatory above cell Peclet 2** -- included so
        // that failure can be measured rather than cited, since it is the
        // scheme the structured Navier-Stokes solvers here still use
        // (DIVIDA_TECNICA.md 4.4).
        Central,
        // Central pulled back towards upwind by a limiter wherever the field
        // is not smooth: second order where that is safe, bounded where it is
        // not. See ConvectionLimiter.hpp.
        LimitedLinearUpwind,
    };

    ImplicitConvectionDiffusionSolver1D(std::size_t nx, double length, double velocity, double diffusivity,
                                         double source, double leftValue, double rightValue);

    // Per-cell coefficient fields; both default to the constant values
    // given to the constructor. Face values are the average of the two
    // adjacent cells (the boundary cell's own value at a domain boundary),
    // matching the edge-averaging idiom the 3D turbulence solvers already
    // use for nu_t. Throws std::invalid_argument on a size mismatch.
    void setVelocityField(const std::vector<double>& velocityPerCell);
    void setDiffusivityField(const std::vector<double>& diffusivityPerCell);

    // Applies to solveBiCGStab() and solveGmres(); ignored by
    // solveGaussSeidel(), which is a stationary method with no Krylov
    // subspace to precondition. Changing it invalidates any cached
    // factorization, which is rebuilt on the next solve.
    void setPreconditioner(Preconditioner preconditioner);

    // Defaults to FirstOrderUpwind, which is what this class has always done
    // and what its existing tests are calibrated against. A non-upwind scheme
    // makes every solve below run outer deferred-correction sweeps; the inner
    // Krylov solve is unchanged.
    void setConvectionScheme(ConvectionScheme scheme) { scheme_ = scheme; }

    // Largest overshoot outside the range spanned by the two Dirichlet
    // values, relative to that range. **Zero is not a tolerance but a
    // theorem**: this equation obeys a maximum principle, so a solution
    // leaving that range is the scheme failing, not the physics. It is the
    // crisp way to catch central differencing above cell Peclet 2, which
    // otherwise only shows up as a suspiciously large error.
    double maxBoundednessViolation() const;

    // |u| h / Gamma at the worst face: the number that decides whether
    // central differencing is admissible at all.
    double maxCellPeclet() const;

    // Reference solver: Gauss-Seidel sweeps directly on the same
    // upwind+central discretization, using each row's true diagonal
    // coefficient (extracted from the operator itself -- see
    // extractTridiagonal()) as the relaxation denominator. Returns
    // iterations actually run.
    std::size_t solveGaussSeidel(std::size_t maxIterations = 20000, double tolerance = 1e-12);

    // BiCGSTAB (Van der Vorst 1992), matrix-free, on the first genuinely
    // non-symmetric operator this project has assembled -- Conjugate
    // Gradient's SPD requirement isn't met here (upwinding makes the two
    // neighbor coefficients unequal). Returns iterations actually run
    // (>=0), or -1 if the algorithm's own breakdown condition (a near-zero
    // denominator) was detected and the solve aborted rather than dividing
    // by ~0.
    long long solveBiCGStab(std::size_t maxIterations = 500, double tolerance = 1e-10);

    // GMRES(m) -- restarted Generalized Minimal RESidual (Saad & Schultz
    // 1986), the other half of the roadmap's long-standing
    // "GMRES/BiCGSTAB" item, on the same non-symmetric operator. Arnoldi
    // with *modified* Gram-Schmidt (deliberately not classical
    // Gram-Schmidt, which loses basis orthogonality catastrophically in
    // floating point), plus Givens rotations applied incrementally to the
    // Hessenberg matrix, so the residual norm is known exactly at every
    // inner step without re-solving the least-squares problem.
    //
    // **What GMRES has that BiCGSTAB does not, and why it earns its extra
    // cost**: GMRES minimizes the residual over the entire Krylov subspace
    // at every step, so its residual norm is monotonically non-increasing
    // *by construction* -- including across restarts, since each cycle
    // starts from the previous cycle's iterate and the zero correction is
    // always inside the new subspace. BiCGSTAB guarantees no such thing and
    // in practice oscillates. That difference is directly measurable rather
    // than merely assertable, so residualHistory() exposes it and this
    // class's tests measure it (the same "find the test that actually
    // distinguishes the two" discipline used for LES-vs-RANS and
    // DES-vs-SST).
    //
    // `restart` is the m in GMRES(m): the Krylov basis is rebuilt from
    // scratch every m inner iterations, bounding storage at m+1 vectors
    // instead of growing without limit. Returns total inner iterations
    // run, or -1 on an Arnoldi/back-substitution breakdown that is not a
    // converged "happy breakdown".
    long long solveGmres(std::size_t restart = 30, std::size_t maxIterations = 500,
                          double tolerance = 1e-10);

    // Relative residual after each iteration of the most recent
    // solveGmres() or solveBiCGStab() call (empty before either has run;
    // solveGaussSeidel() does not populate it). Exposed specifically so the
    // monotonicity contrast described above can be measured instead of
    // assumed. **With a preconditioner set this is the preconditioned
    // residual ||M^-1(b - A x)|| / ||M^-1 b||**, which is what those
    // solvers actually minimize -- not the same quantity as the
    // unpreconditioned residual, so histories from different preconditioner
    // settings are not directly comparable in magnitude.
    const std::vector<double>& residualHistory() const { return residualHistory_; }

    double value(std::size_t i) const { return phi_[i]; }
    double cellCenterX(std::size_t i) const;

private:
    // Ghost-mirrored Dirichlet lookup: x[i] for 0<=i<nx, the mirrored wall
    // value for i==-1 or i==nx.
    double phiAt(const std::vector<double>& x, long long i) const;

    // Face-centered coefficients; face f sits between cells f-1 and f, so
    // f runs 0..nx (0 = left wall, nx = right wall).
    double velocityAtFace(std::size_t face) const;
    double diffusivityAtFace(std::size_t face) const;

    // The full (affine) finite-volume residual at cell i: (convective flux
    // out - in) - (diffusive flux out - in), using the true ghost values
    // (which carry the Dirichlet boundary constants). Solving
    // stencilAt(phi, i) == source*h for every i is exactly this class's
    // linear system.
    double stencilAt(const std::vector<double>& x, std::size_t i) const;

    // Convective face value under the active scheme, and the per-cell
    // difference between that and the upwind value the matrix assembles --
    // the deferred correction, ready to be added to the right-hand side.
    double convectedFaceValue(const std::vector<double>& x, std::size_t face) const;
    std::vector<double> convectionCorrection(const std::vector<double>& x) const;

    // Runs `innerSolve` repeatedly, refreshing the deferred correction from
    // the previous iterate between sweeps, until the solution stops moving.
    // With FirstOrderUpwind there is no correction, so it runs exactly once
    // and this class behaves exactly as it always did.
    void outerSweeps(const std::function<void()>& innerSolve, double tolerance);

    // The single-sweep bodies; the public methods above wrap these in the
    // deferred-correction loop.
    std::size_t solveGaussSeidelOnce(std::size_t maxIterations, double tolerance);
    long long solveBiCGStabOnce(std::size_t maxIterations, double tolerance);
    long long solveGmresOnce(std::size_t restart, std::size_t maxIterations, double tolerance);

    // The pure linear part of stencilAt (the boundary constants subtracted
    // out) -- what a Krylov method's matrix-vector product must be, since
    // Krylov search directions have no meaning for an affine map's constant
    // term.
    std::vector<double> applyOperator(const std::vector<double>& x) const;

    // source*h minus the boundary constants stencilAt's ghosts contribute
    // at the two boundary cells -- so that applyOperator(phi) == rhs() is
    // exactly the same system stencilAt(phi,i) == source*h describes.
    std::vector<double> rightHandSide() const;

    // Recovers the three matrix diagonals by probing applyOperator() with
    // unit vectors, one column at a time -- deliberately NOT by
    // hand-transcribing the coefficients out of stencilAt()'s formula.
    // Two reasons: (1) the boundary rows differ from the interior ones
    // (the ghost mirror folds a neighbor coefficient back onto the
    // diagonal with the opposite sign), which is exactly the sort of
    // detail a hand-derived duplicate gets subtly wrong; (2) this project
    // has already been bitten by precisely that failure mode --
    // MixingLengthChannelFlowSolver1D's frictionVelocity() once assumed a
    // wall treatment its own solve() did not actually use, a ~38% error.
    // Probing costs nx matvecs once per solve, which is setup-time noise
    // here, and it is correct by construction rather than by review.
    //
    // Also *verifies* tridiagonality while extracting (every probed column
    // must be zero outside rows j-1, j, j+1), turning the assumption the
    // ILU(0) reasoning rests on into a checked fact.
    void extractTridiagonal();

    // Builds the ILU(0) factors from the extracted diagonals. For a
    // tridiagonal matrix this is the exact (Thomas-algorithm) LU
    // factorization -- see the Preconditioner enum's comment.
    void buildIncompleteLU();

    // M^-1 * v for the currently selected preconditioner.
    std::vector<double> applyPreconditioner(const std::vector<double>& v) const;

    // Prepares whatever the current preconditioner needs, once per solve.
    void prepareSolve();

    static double dot(const std::vector<double>& a, const std::vector<double>& b);

    std::size_t nx_;
    double length_;
    double h_;
    std::vector<double> velocity_;    // per cell
    std::vector<double> diffusivity_; // per cell
    double source_;
    double leftValue_;
    double rightValue_;
    std::vector<double> phi_;
    std::vector<double> residualHistory_;

    Preconditioner preconditioner_ = Preconditioner::None;
    ConvectionScheme scheme_ = ConvectionScheme::FirstOrderUpwind;
    // Refreshed between outer sweeps; empty (and unused) under upwind.
    std::vector<double> convectionCorrection_;
    // Extracted matrix diagonals: sub_[i] multiplies phi[i-1], diag_[i]
    // multiplies phi[i], super_[i] multiplies phi[i+1].
    std::vector<double> sub_;
    std::vector<double> diag_;
    std::vector<double> super_;
    // ILU(0) factors: iluDiag_ is U's diagonal, iluSub_ is L's subdiagonal
    // (U's superdiagonal is super_ unchanged, L's diagonal is 1).
    std::vector<double> iluDiag_;
    std::vector<double> iluSub_;
};

} // namespace aether::solver
