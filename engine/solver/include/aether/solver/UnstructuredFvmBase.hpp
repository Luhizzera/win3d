#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

namespace aether::solver {

// Shared finite-volume machinery for every solver in this project that runs
// on a tetrahedral mesh, extracted from UnstructuredDiffusionSolver and
// UnstructuredCavitySolver3D.
//
// **Why it was extracted.** The two classes had grown independent copies of
// the same five things: the over-relaxed face decomposition, the face
// coefficient |A|^2/(A.d), the inverse-distance least-squares gradient
// stencil, the symmetric 3x3 inverse with its rank guard, and the
// matrix-free Conjugate Gradient loop. That is the same situation that
// forced StaggeredCavityBase3D out of the six 3D cavities, and the reason
// recorded there applies unchanged: **porting a correction N times is the
// expensive version of the problem.** It was not a hypothetical risk here --
// the copies had already diverged within a day of being written, in the
// three ways recorded as DIVIDA_TECNICA.md 1.2, 1.3 and 1.4.
//
// **All three of those divergences are now resolved here**, which is what
// the extraction was for: the Laplacian, its non-orthogonal deferred
// correction, and the gradient (least-squares with a Green-Gauss fallback,
// never a silent zero) exist once and both solvers use them. Where the two
// solvers genuinely differ, the difference is a named argument at the call
// site -- which boundary faces enter the operator, what value a boundary
// face carries, whether a reference cell is pinned -- rather than a silent
// divergence between two bodies of code nobody diffs.
//
// **No virtual functions**, for the same reason StaggeredCavityBase3D has
// none: the customisation points are build-time (std::function, called once
// per face while assembling, or stored and consulted per boundary face) and
// the hot path is a template (conjugateGradient), so nothing in an inner
// loop pays for dispatch. The destructor is protected and non-virtual --
// this is a base to reuse code, not to delete through.
class UnstructuredFvmBase {
public:
    std::size_t cellCount() const { return mesh_->cellCount(); }

    // Largest |A_nonorth| / |A| over the interior faces: 0 on a perfectly
    // orthogonal mesh, approaching 1 as a face becomes parallel to the line
    // between the cell centres. Exposed because it is the quantity that
    // decides whether the non-orthogonal correction matters, and a caller
    // measuring accuracy should be able to report it alongside.
    double maxNonOrthogonality() const;

    // How many cells have a least-squares stencil too rank-deficient to pin
    // all three gradient components -- a mesh-quality diagnostic in the same
    // sense as maxNonOrthogonality(), and the one that says how much work
    // the Green-Gauss fallback is doing. Zero means the fallback never runs
    // on this mesh; a large count means the gradients are largely
    // first-order and the solver's formal accuracy claim does not hold.
    std::size_t deficientStencilCount() const { return deficientStencilCount_; }

    // How a convected quantity is reconstructed at a face.
    enum class ConvectionScheme {
        // phi_f = phi_upwind. Bounded unconditionally and first-order
        // accurate. Its error is numerical diffusion, and the part that makes
        // it more than a nuisance is that the diffusion is largest for flow
        // *oblique* to the face -- which on tetrahedra is every face.
        FirstOrderUpwind,
        // Second-order reconstruction from the upwind cell, pulled back
        // towards upwind by a limiter wherever the solution is not smooth.
        // See limitedFaceValue().
        LimitedLinearUpwind,
    };

protected:
    explicit UnstructuredFvmBase(const mesh::TetrahedralMesh& mesh) : mesh_(&mesh) {}
    ~UnstructuredFvmBase() = default;

    // -- Face geometry ---------------------------------------------------

    // The over-relaxed decomposition of the face area vector A about the
    // centroid-to-centroid vector d:
    //
    //   A_orth = (|A|^2 / (A . d)) d,      A_nonorth = A - A_orth
    //
    // so a diffusive flux through the face splits into an implicit part
    // a_f (phi_N - phi_P), with a_f = |A|^2/(A.d) > 0, plus a correction
    // (grad phi)_f . A_nonorth that no two-point stencil can represent.
    // Both solvers need exactly this, which is why it lives here.
    struct InteriorFace {
        std::size_t owner;
        std::size_t neighbour;
        double coefficient;              // a_f = |A|^2 / (A . d)
        core::Vector3 areaVector;        // owner -> neighbour
        core::Vector3 nonOrthogonalArea; // A - a_f d, as seen from the owner
        core::Vector3 unitD;             // normalised c_N - c_P
        double distance;                 // |c_N - c_P|
        // Interpolation weight for the owner, from the perpendicular
        // distances of the two centroids to the face plane rather than a
        // flat 0.5: on a skewed mesh the two cells sit at different
        // distances, and 0.5 there is a first-order error in every face
        // quantity built from it.
        double ownerWeight;
    };

    // A boundary face's "neighbour" is the face centroid itself, so d_b runs
    // from the cell centroid to the face centroid and the same decomposition
    // applies. Every boundary face gets an entry, whatever its condition:
    // what the condition changes is only whether the coefficient enters the
    // operator.
    struct BoundaryFace {
        std::size_t meshFace; // index into the mesh's face list
        std::size_t cell;
        double coefficient;              // a_b = |A|^2 / (A . d_b)
        core::Vector3 areaVector;        // outward
        core::Vector3 centroid;
        core::Vector3 nonOrthogonalArea;
        core::Vector3 unitD;
        double distance;
        bool entersOperator = false; // see buildFaceGeometry
    };

    // Fills interiorFaces_, boundaryFaces_ and the two diagonals in a single
    // pass over the mesh.
    //
    // `boundaryEntersOperator(meshFace)` decides whether a boundary face
    // contributes its coefficient to the Laplacian's diagonal. That is where
    // the two solvers' assembly genuinely differs, and both readings are the
    // same rule stated for a different condition: a **prescribed value**
    // (Dirichlet -- a fixed temperature, an outlet pressure) makes the face an
    // extra connection and belongs in the operator, while a **zero-gradient**
    // face (an insulated wall, the pressure condition on a solid wall)
    // carries no flux and therefore no coefficient at all. Passing {} means
    // no boundary face enters, which leaves the operator singular unless a
    // reference cell is pinned.
    //
    // Faces with A . d <= 0 are skipped entirely rather than assembled with a
    // negative coefficient: that would mean the outward normal points back
    // towards the cell centre, i.e. a broken mesh, and a negative coefficient
    // destroys positive definiteness silently -- the CG below would then
    // converge to something plausible and wrong.
    void buildFaceGeometry(const std::function<bool(std::size_t)>& boundaryEntersOperator);

    // -- Boundary values ---------------------------------------------------

    // Registers what a boundary face's value *is*, when it has one:
    // `boundaryFaceValue(meshFace, value)` returns true and writes `value`
    // for a face carrying a prescribed value, false for a zero-gradient one.
    //
    // Two different questions are answered from this one callback, and both
    // need it: Green-Gauss needs a value on every face of the closed surface
    // integral (a zero-gradient face contributes the cell's own value), and
    // the non-orthogonal correction needs it on the faces that entered the
    // operator. Keeping it as one registered function is what stops the two
    // solvers from answering the same question differently again.
    //
    // **Evaluated once per boundary face and cached, never called from
    // inside a solve.** It used to be called from the Green-Gauss loop and
    // from the non-orthogonal correction, i.e. thousands of times per step,
    // which contradicted this class's own claim that no inner loop pays for
    // dispatch. Writing the Python bindings is what made that visible and
    // not merely inelegant: a predicate defined in Python would have taken
    // the GIL on every one of those calls, which would have made the bindings
    // too slow for the experiments they exist to enable.
    void setBoundaryFaceValue(std::function<bool(std::size_t, double&)> boundaryFaceValue) {
        boundaryFaceValue_ = std::move(boundaryFaceValue);
        refreshBoundaryValues();
    }

    // Re-evaluates the registered callback for every boundary face. Called
    // automatically by buildFaceGeometry() and by setBoundaryFaceValue(), so
    // the two possible orderings both end up with a filled cache; a solver
    // whose prescribed values change without a rebuild must call it.
    void refreshBoundaryValues();

    // -- Gradients ---------------------------------------------------------

    static constexpr std::size_t kBoundaryStencil = static_cast<std::size_t>(-1);

    struct GradientStencilEntry {
        std::size_t neighbour; // kBoundaryStencil when the "neighbour" is a prescribed face value
        double boundaryValue;
        core::Vector3 weightedDelta; // w_i * d_i
    };

    // Inverse of the symmetric 3x3 normal-equation matrix, stored as its six
    // unique components. Geometry only, so it is built once alongside the
    // face coefficients.
    struct SymmetricInverse {
        double xx = 0.0, xy = 0.0, xz = 0.0, yy = 0.0, yz = 0.0, zz = 0.0;
        bool valid = false; // false => the stencil was rank-deficient
        core::Vector3 apply(const core::Vector3& v) const {
            return {xx * v.x + xy * v.y + xz * v.z, xy * v.x + yy * v.y + yz * v.z,
                    xz * v.x + yz * v.y + zz * v.z};
        }
    };

    // Inverse of a symmetric 3x3 matrix given by its six unique components,
    // via cofactor expansion. `invalid` (false) when the matrix is singular
    // relative to its own magnitude -- scaled by scale^3 so the test is
    // dimensionless rather than an absolute threshold that would depend on
    // the mesh's units. Shared by every symmetric 3x3 solve in this class:
    // the least-squares gradient stencil below, and the flux-reconstruction
    // stencil UnstructuredCavitySolver3D builds on top of this base.
    static SymmetricInverse invertSymmetric3x3(double xx, double xy, double xz, double yy, double yz,
                                                double zz);

    // Builds the inverse-distance-weighted least-squares stencil: the fit
    // that best reproduces the measured differences to every neighbour,
    // minimising sum_i w_i (grad(phi)_P . d_i - (phi_i - phi_P))^2 with
    // w_i = 1/|d_i|^2. Unlike Green-Gauss it stays second-order on a skewed
    // mesh, which is why it exists: the Green-Gauss version's convergence
    // order stalled around 0.6-1.0 instead of reaching 2 (Fase 2.2).
    //
    // `useBoundaryValues` includes the registered boundary values as extra
    // stencil entries. A face with no prescribed value is never included
    // either way: a zero-gradient face asserts that the *normal derivative*
    // vanishes, which this fit would misread as "the value there equals the
    // cell value", biasing the tangential components too. Requires
    // setBoundaryFaceValue() to have been called when true.
    void buildGradientStencils(bool useBoundaryValues);

    // The gradient this project uses: least squares where the stencil
    // supports it, Green-Gauss where it does not.
    //
    // **The fallback is the point, and it is not optional.** One of the two
    // copies this base replaced returned a *zero* gradient for a
    // rank-deficient cell, which is not "no answer" but a wrong answer that
    // looks plausible -- and deficient stencils appear exactly where the
    // geometry is complicated and the gradient matters most. That was
    // DIVIDA_TECNICA.md 1.3; there is now no way for a caller to make that
    // choice, because the choice is not offered. deficientStencilCount()
    // says how often the fallback runs, so its cost to accuracy stays
    // visible rather than silent.
    // `clampFallback` gates the steepest-slope bound on deficient-stencil
    // cells (see the .cpp for why the bound exists). Default true for every
    // ordinary caller. A caller building a genuinely *linear* operator out
    // of this function -- the GMRES coupling correction in
    // UnstructuredCavitySolver3D is the only one so far -- must pass false:
    // the bound is a min() against a slope computed from `field` itself, so
    // with it left on this function is not linear, not even homogeneous of
    // degree 1, and a Krylov method silently converges to the wrong answer
    // if handed an operator that lies about being linear.
    std::vector<core::Vector3> computeCellGradients(const std::vector<double>& field,
                                                     bool clampFallback = true) const;

    // Green-Gauss: grad(phi)_P = (1/V_P) sum_f phi_f A_f. Cheap, but only
    // first-order on a skewed mesh, which is why it is the fallback and not
    // the default.
    std::vector<core::Vector3> greenGaussGradients(const std::vector<double>& field) const;

    // Distance-weighted interpolation of cell gradients to a face. Both the
    // non-orthogonal correction and the Rhie-Chow face flux start from this
    // and then replace its component along d by the compact difference
    // (phi_N - phi_P)/|d| -- the one direction two cell values pin down
    // accurately. They apply that replacement with opposite signs, so the
    // shared part stops here rather than being forced into one expression.
    core::Vector3 averagedFaceGradient(const InteriorFace& face,
                                        const std::vector<core::Vector3>& gradients) const {
        return gradients[face.owner] * face.ownerWeight +
               gradients[face.neighbour] * (1.0 - face.ownerWeight);
    }

    // -- Convection ---------------------------------------------------------

    // Face value of `field` for a convective flux across `face`, given the
    // mass flux through it and the cell gradients.
    //
    // **The limited form, and why this shape of it.** Writing the face value
    // as a correction to the upwind value,
    //
    //   phi_f = phi_C + psi * (phi_central - phi_C)
    //
    // makes psi = 0 exactly upwind and psi = 1 exactly central, so the
    // limiter interpolates between a scheme that is bounded and a scheme that
    // is second order. `phi_central` is the distance-weighted interpolation
    // the rest of this class already uses, not a flat average -- on a skewed
    // mesh the two differ at first order, which would defeat the point.
    //
    // The limiter argument is the standard unstructured ratio
    //
    //   r = 2 (grad(phi)_C . d_CD) / (phi_D - phi_C) - 1
    //
    // which compares the *upwind cell's own gradient* against the difference
    // it is being asked to span: r = 1 in a smooth field (psi = 1, second
    // order), r <= 0 at an extremum (psi = 0, upwind, bounded). The gradient
    // is the least-squares one this class already computes, which is what
    // makes second-order-accurate reconstruction possible without a wider
    // stencil -- there is no "upwind-upwind" cell to reach for on tetrahedra.
    //
    // van Leer's psi(r) = (r + |r|) / (1 + |r|): smooth, symmetric, and
    // inside the TVD region. Smoothness matters here beyond taste -- a
    // limiter with a kink (minmod, superbee) makes the residual of a steady
    // iteration non-differentiable and can stall it short of convergence.
    double faceValue(const InteriorFace& face, double massFlux, const std::vector<double>& field,
                      const std::vector<core::Vector3>& gradients, ConvectionScheme scheme) const;

    // -- The Laplacian -----------------------------------------------------

    static constexpr std::size_t kNoPinnedCell = static_cast<std::size_t>(-1);

    // A x, where row P reads
    //   (sum of all face coefficients touching P) x_P  -  sum_f a_f x_N
    // Symmetric by construction: each interior face contributes the same a_f
    // to both (P,N) and (N,P), which is what lets plain CG solve it.
    //
    // `pinnedCell` replaces one row by the identity, removing the null space
    // of a pure-Neumann problem (every boundary a wall). Pass kNoPinnedCell
    // when a Dirichlet face already fixes the level -- pinning as well would
    // over-constrain the system with two incompatible references.
    std::vector<double> applyLaplacian(const std::vector<double>& x,
                                        std::size_t pinnedCell = kNoPinnedCell) const;

    // The non-orthogonal flux the implicit part cannot represent, evaluated
    // from the current iterate and accumulated per cell, ready to be added to
    // the right-hand side.
    //
    // **Not optional, and that was measured rather than assumed.** Fase 2.2
    // ran the mesh convergence study with only the implicit part: rms error
    // 2.879 -> 2.433 -> 2.365 for n = 4/6/8, i.e. observed order 0.42 then
    // 0.10 -- a plateau, not convergence. Non-orthogonality is a property of
    // cell *shape*, so refining a tet mesh does not reduce it (~1.5 on these
    // meshes) and the dropped term settles to a mesh-quality error floor
    // instead of vanishing.
    std::vector<double> nonOrthogonalCorrection(const std::vector<double>& field) const;

    // -- Linear solve ------------------------------------------------------

    static double dot(const std::vector<double>& a, const std::vector<double>& b);

    // Solves the full Laplacian -- implicit part plus non-orthogonal
    // correction -- by **deferred correction**: each outer sweep evaluates
    // the correction from the previous iterate, moves it to the right-hand
    // side, and runs a plain CG on the implicit part. The matrix therefore
    // stays a symmetric positive-definite M-matrix and nothing about the CG
    // has to change, which is the whole reason for doing it this way.
    //
    // `x` is the initial guess and the result. Returns the number of outer
    // sweeps run and writes the largest per-cell change across the final one
    // to `lastOuterChange` -- the number that distinguishes a solve that
    // settled from one that merely hit the sweep cap, which mean very
    // different things about the result.
    // `relaxation` is in/out and starts at 1.0. It is a parameter rather than
    // a local because a solver that runs this every time step needs it to
    // **persist**: the damping is learned from watching the iteration
    // misbehave, and a projection that allows only four sweeps per step never
    // gets far enough to learn it if the factor resets each time. Carrying it
    // means the mesh is diagnosed once and the answer kept.
    std::size_t solveDeferredCorrection(std::vector<double>& x, const std::vector<double>& baseRhs,
                                         std::size_t pinnedCell, std::size_t maxIterations,
                                         double tolerance, std::size_t maxOuterSweeps,
                                         double& lastOuterChange, double& relaxation) const;

    // Matrix-free Conjugate Gradient. `apply(v)` returns A v; `x` is both the
    // initial guess and the result, which matters -- the pressure Poisson
    // solve starts from the previous step's pressure and converges in a
    // handful of iterations because of it.
    //
    // A template rather than a std::function because this is the one inner
    // loop here: the operator is applied once per iteration over every cell,
    // and an indirect call on that path would be paid for no benefit. `A`
    // must be symmetric positive definite -- everything assembled by
    // buildFaceGeometry is, by construction, provided at least one boundary
    // face enters the operator or a reference cell is pinned.
    template <typename ApplyOperator>
    static void conjugateGradient(ApplyOperator&& apply, const std::vector<double>& rhs,
                                   std::vector<double>& x, std::size_t maxIterations, double tolerance) {
        const std::size_t n = rhs.size();
        std::vector<double> residual(n);
        {
            const std::vector<double> applied = apply(x);
            for (std::size_t i = 0; i < n; ++i) {
                residual[i] = rhs[i] - applied[i];
            }
        }
        std::vector<double> direction = residual;
        double residualDotResidual = dot(residual, residual);

        for (std::size_t iteration = 0; iteration < maxIterations; ++iteration) {
            if (std::sqrt(residualDotResidual) < tolerance) {
                break;
            }
            const std::vector<double> applied = apply(direction);
            const double directionDotApplied = dot(direction, applied);
            if (directionDotApplied == 0.0) {
                break;
            }
            const double alpha = residualDotResidual / directionDotApplied;
            double newResidualDotResidual = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                x[i] += alpha * direction[i];
                residual[i] -= alpha * applied[i];
                newResidualDotResidual += residual[i] * residual[i];
            }
            const double beta = newResidualDotResidual / residualDotResidual;
            for (std::size_t i = 0; i < n; ++i) {
                direction[i] = residual[i] + beta * direction[i];
            }
            residualDotResidual = newResidualDotResidual;
        }
    }

    // **Restarted GMRES(m), matrix-free, for operators conjugateGradient()
    // cannot solve.** Extracted from ImplicitConvectionDiffusionSolver1D's
    // own (independently tested) implementation rather than rewritten --
    // this project's own lesson from item 2.1 is that porting a correction N
    // times is the expensive version of the problem, and re-deriving a
    // Krylov method's numerics from scratch is exactly the kind of thing
    // worth not doing twice. The 1D solver keeps its own copy for now rather
    // than being refactored onto this one in the same change that adds it;
    // that duplication is real and is recorded, not silently left
    // (DIVIDA_TECNICA.md 4.3's follow-up note).
    //
    // **Why this exists at all**, when the base already has
    // conjugateGradient(): CG requires `A` symmetric positive definite.
    // UnstructuredCavitySolver3D's pressure-velocity coupling on a distorted
    // mesh is not -- DIVIDA_TECNICA.md 4.3 measured a scalar-relaxed defect
    // correction failing on exactly this operator, including on meshes where
    // the baseline algorithm already worked, because the operator's
    // contraction varies by *direction* in residual space and a single
    // relaxation number cannot describe direction-dependent behaviour.
    // Modified Gram-Schmidt Arnoldi plus incremental Givens rotations search
    // the whole Krylov subspace at each step instead of scaling one guessed
    // direction, which is what a genuinely non-symmetric operator needs.
    //
    // `x` is both the initial guess and the result. Returns the number of
    // matrix-vector products used. Like conjugateGradient(), this does not
    // itself signal whether the result actually converged -- a caller that
    // needs to know checks the residual, the same way it would for CG; both
    // stay simple matrix-free primitives rather than carrying a policy for
    // what "not converged" should mean to every different caller.
    template <typename ApplyOperator>
    static std::size_t gmres(ApplyOperator&& apply, const std::vector<double>& rhs,
                              std::vector<double>& x, std::size_t restart, std::size_t maxIterations,
                              double tolerance) {
        const std::size_t n = rhs.size();
        const double bNorm = std::max(std::sqrt(dot(rhs, rhs)), 1e-300);
        const std::size_t m = std::max<std::size_t>(restart, 1);
        const double epsilon = std::numeric_limits<double>::epsilon();
        std::size_t totalIterations = 0;

        while (totalIterations < maxIterations) {
            // --- Each restart cycle begins afresh from the current iterate,
            // which is why restarted GMRES stays monotone globally: the zero
            // correction is always available inside the new Krylov
            // subspace, so a cycle can never make the residual worse than it
            // started.
            std::vector<double> r(n);
            {
                const std::vector<double> ax = apply(x);
                for (std::size_t i = 0; i < n; ++i) {
                    r[i] = rhs[i] - ax[i];
                }
            }
            const double beta = std::sqrt(dot(r, r));
            if (beta / bNorm < tolerance) {
                return totalIterations;
            }

            std::vector<std::vector<double>> basis; // Arnoldi basis vectors v_0..v_k
            basis.reserve(m + 1);
            basis.push_back(r);
            for (double& value : basis[0]) {
                value /= beta;
            }

            // hessenberg[j] is column j of the (rotated) Hessenberg matrix,
            // holding its j+2 potentially-nonzero entries.
            std::vector<std::vector<double>> hessenberg;
            std::vector<double> cosines(m, 0.0);
            std::vector<double> sines(m, 0.0);
            std::vector<double> g(m + 1, 0.0); // rotated RHS of the least-squares problem
            g[0] = beta;

            std::size_t completed = 0;
            bool converged = false;
            bool breakdown = false;

            for (std::size_t j = 0; j < m && totalIterations < maxIterations; ++j) {
                std::vector<double> w = apply(basis[j]);

                // Modified Gram-Schmidt: each projection is subtracted from
                // the *already-updated* w rather than the original --
                // mathematically identical to classical Gram-Schmidt,
                // numerically far better, the standard choice inside Arnoldi
                // for exactly this reason.
                std::vector<double> column(j + 2, 0.0);
                for (std::size_t i = 0; i <= j; ++i) {
                    column[i] = dot(w, basis[i]);
                    for (std::size_t k = 0; k < n; ++k) {
                        w[k] -= column[i] * basis[i][k];
                    }
                }
                const double arnoldiNorm = std::sqrt(dot(w, w));
                column[j + 1] = arnoldiNorm;

                // Replay every previous rotation onto the new column, then
                // build the one that annihilates its subdiagonal entry --
                // keeps an up-to-date QR factorization of the Hessenberg
                // matrix incrementally, so |g[j+1]| below *is* the exact
                // residual norm of the least-squares solution, available
                // without solving it.
                for (std::size_t i = 0; i < j; ++i) {
                    const double rotated = cosines[i] * column[i] + sines[i] * column[i + 1];
                    column[i + 1] = -sines[i] * column[i] + cosines[i] * column[i + 1];
                    column[i] = rotated;
                }
                const double denominator = std::hypot(column[j], column[j + 1]);
                if (denominator <= epsilon * beta) {
                    breakdown = true;
                    break;
                }
                cosines[j] = column[j] / denominator;
                sines[j] = column[j + 1] / denominator;
                column[j] = denominator;
                column[j + 1] = 0.0;
                g[j + 1] = -sines[j] * g[j];
                g[j] = cosines[j] * g[j];

                hessenberg.push_back(std::move(column));
                ++completed;
                ++totalIterations;

                const double relativeResidual = std::fabs(g[j + 1]) / bNorm;
                // A zero Arnoldi norm is a "happy breakdown": the Krylov
                // subspace is already invariant, so the least-squares
                // solution over it is the exact solution -- success, not
                // failure.
                if (relativeResidual < tolerance || arnoldiNorm <= epsilon * beta) {
                    converged = true;
                    break;
                }
                for (double& value : w) {
                    value /= arnoldiNorm;
                }
                basis.push_back(std::move(w));
            }

            // Back-substitute the (upper-triangular, post-rotation) system
            // H*y = g, then apply the correction x += V*y. Done at the end
            // of every cycle, converged or not, so a restart resumes from
            // the best iterate this cycle found rather than discarding its
            // work.
            if (completed > 0) {
                std::vector<double> y(completed, 0.0);
                bool solvable = true;
                for (std::size_t row = completed; row-- > 0;) {
                    double sum = g[row];
                    for (std::size_t col = row + 1; col < completed; ++col) {
                        sum -= hessenberg[col][row] * y[col];
                    }
                    if (std::fabs(hessenberg[row][row]) <= epsilon * beta) {
                        solvable = false;
                        break;
                    }
                    y[row] = sum / hessenberg[row][row];
                }
                if (!solvable) {
                    return totalIterations; // keep the best x found so far
                }
                for (std::size_t col = 0; col < completed; ++col) {
                    for (std::size_t k = 0; k < n; ++k) {
                        x[k] += y[col] * basis[col][k];
                    }
                }
            }

            if (breakdown || converged) {
                return totalIterations;
            }
        }
        return totalIterations;
    }

    const mesh::TetrahedralMesh* mesh_;
    std::vector<InteriorFace> interiorFaces_;
    std::vector<BoundaryFace> boundaryFaces_;
    // sum_f a_f per cell over the **interior** faces only. Kept apart from
    // laplacianDiagonal_ because a viscous operator and a pressure Poisson
    // operator do not agree on which boundary faces belong in the diagonal:
    // an outlet is Dirichlet for pressure and zero-gradient for velocity, and
    // sharing one diagonal between the two silently imposed u = 0 there
    // (DIVIDA_TECNICA.md 1.4).
    std::vector<double> interiorDiagonal_;
    // interiorDiagonal_ plus the boundary faces that entered the operator:
    // the diagonal of the Laplacian applyLaplacian() assembles.
    std::vector<double> laplacianDiagonal_;
    std::vector<std::vector<GradientStencilEntry>> gradientStencil_;
    std::vector<SymmetricInverse> gradientMatrixInverse_;
    std::function<bool(std::size_t, double&)> boundaryFaceValue_;
    // The callback's answer, cached per *mesh* face -- mesh-indexed rather
    // than boundary-indexed because Green-Gauss integrates over a cell's
    // whole closed surface, mesh face by mesh face, including any face
    // buildFaceGeometry() rejected. Dropping a face there would break the
    // divergence-theorem identity the method rests on.
    std::vector<char> boundaryHasValue_;
    std::vector<double> boundaryValueCache_;
    std::size_t deficientStencilCount_ = 0;
};

} // namespace aether::solver
