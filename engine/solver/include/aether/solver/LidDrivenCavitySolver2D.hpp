#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Solves the 2D incompressible Navier-Stokes equations on a square cavity
// with solid walls: no-slip (u=v=0) on the left, right and bottom walls,
// and a moving lid on top -- the classic lid-driven cavity benchmark. Same
// numerical method as TaylorGreenVortexSolver2D (collocated grid, Chorin
// projection: explicit convection+diffusion predictor, then a Conjugate-
// Gradient pressure projection), but with solid-wall ghost-cell boundary
// treatment instead of periodic wraparound: velocity ghost values are
// mirrored to enforce the exact wall value (Dirichlet, e.g.
// ghost = 2*wallValue - interior), pressure ghost values simply copy the
// interior value (Neumann/zero-gradient, the standard assumption at an
// impermeable wall with no-slip velocity).
//
// The lid's tangential velocity is, by default, *regularized* to taper
// smoothly to zero at its two top corners (lidVelocity * sin^2(pi*x/lengthX))
// rather than holding a constant value with a step discontinuity where it
// meets the stationary side walls. A discontinuous lid is a genuinely
// singular boundary condition for the continuous Navier-Stokes equations
// (pressure blows up at those two corners); measured directly while
// validating this class, that singularity dominated the discrete divergence
// (an interior cell one step away from a top corner had ~20x the divergence
// of any other interior cell) rather than it being spread out as ordinary
// discretization error. Tapering the lid is standard, documented practice
// for exactly this reason, not an ad hoc fix.
//
// **`taperLid = false` gives the classical uniform lid instead** (constant
// lidVelocity right up to both corners), added specifically to let a caller
// reproduce Ghia, Ghia & Shin (1982)'s own posed problem rather than this
// class's regularized one -- see
// python/research/ghia_1982_validation.py, which measured the two as
// genuinely different problems (refining the mesh did not close the gap
// against the tapered lid's results) before this option existed to test
// why directly. Defaults to the taper: every number this class has ever
// produced was measured with it, and every existing test of this class
// keeps passing unchanged with the default untouched.
//
// This class was originally validated *without* comparing to literature
// benchmark tables (e.g. Ghia et al. 1982): recalling one from memory would
// have been a real accuracy risk. It has since been compared for real (see
// the file referenced above), with the reference data fetched from actual
// sources rather than memory -- validated instead, at the time, against
// properties derivable from first principles: the exact rest state when the
// lid isn't moving, mass conservation (bounded discrete divergence, as with
// TaylorGreenVortexSolver2D), and the flow-reversal that mass conservation
// in a closed box necessarily forces (see solver_tests.cpp for details).
//
// Same collocated-grid limitation as TaylorGreenVortexSolver2D applies:
// checkerboard-prone without Rhie-Chow interpolation (not implemented).
class LidDrivenCavitySolver2D {
public:
    // How the convected velocity is reconstructed at a face.
    //
    // **The default is the limited scheme, and it was measured before being
    // made the default** -- the whole content of DIVIDA_TECNICA.md 4.4.
    //
    // Central differencing is unbounded above cell Reynolds 2, which the 1D
    // case with a known exact solution shows costing 3.17 times the boundary
    // range in overshoot at cell Peclet 8.3
    // (testConvectionSchemesAgainstCellPeclet). Whether that harm actually
    // reaches this solver was a separate question, and the answer needed a
    // reference: the cavity at Re=400 on n=128 (cell Reynolds 3.1), against
    // which the three schemes were run at n=16 and n=32. RMS deviation of u
    // along the vertical centreline:
    //
    //         n=16      n=32     cell Reynolds 25.0 and 12.5
    //   central    0.066457  0.023450
    //   upwind     0.072535  0.040897
    //   limited    0.054422  0.021399    <- best at both
    //
    // So the limited scheme is not merely safer, it is more accurate here --
    // 18% better than central at n=16 and 9% at n=32, and 25% to 48% better
    // than plain upwind. That, plus boundedness central cannot offer, is what
    // justified changing a validated solver's default.
    //
    // Central is kept because every number this solver produced before this
    // change was measured with it, and reproducing them has to stay possible.
    // The six 3D closures built on StaggeredCavityBase3D are **still
    // central**: they share a different predictor, and porting it is its own
    // measurement.
    enum class ConvectionScheme {
        // The original: u . grad(u) with central differences, non-conservative.
        Central,
        // Conservative face fluxes with the upwind face value. Bounded at any
        // cell Reynolds number, first-order.
        FirstOrderUpwind,
        // Conservative face fluxes with the limited linear-upwind value --
        // second order where the field is smooth, upwind where it is not. Uses
        // the same limiter as the unstructured side (ConvectionLimiter.hpp).
        LimitedLinearUpwind,
    };

    LidDrivenCavitySolver2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                             double viscosity, double lidVelocity,
                             ConvectionScheme convection = ConvectionScheme::LimitedLinearUpwind,
                             bool taperLid = true);

    // The explicit stability limit **with the safety factor applied**, so
    // the returned value may be stepped with directly -- see
    // ExplicitTimeStep.hpp for why that is worth stating. It used to return
    // the marginal value, with every caller writing its own `0.3 *` in front
    // and two other engine layers not writing one at all
    // (DIVIDA_TECNICA.md 4.1).
    //
    // The margin this buys is against *perturbation*, which is what Fase 1
    // measured biting: a modest change to the scheme took this case to NaN in
    // ~500 steps at CFL 1.0000. It is not a fix for the other half of that
    // item -- central differencing of convection at cell Reynolds 16.7, which
    // is a property of the scheme and not of the step, and which a smaller dt
    // cannot cure.
    double stableTimeStep() const;

    void step(double dt);

    double u(std::size_t i, std::size_t j) const;
    double v(std::size_t i, std::size_t j) const;
    double pressure(std::size_t i, std::size_t j) const;
    double time() const { return time_; }

    // Overwrites the full state (velocity, pressure, simulated time) with
    // externally supplied data -- the counterpart to reading u()/v()/
    // pressure()/time() out to save a checkpoint (see
    // engine/persistence/FieldArchive). Throws std::invalid_argument if any
    // field's size doesn't match nx*ny. Deliberately not validated any
    // further than size (e.g. no divergence check): the caller is
    // responsible for supplying a state that actually came from a solver
    // with the same grid, the same way a constructor doesn't validate the
    // physical plausibility of its parameters either.
    void loadState(std::vector<double> u, std::vector<double> v, std::vector<double> p, double time);

    double maxDivergence() const;

    // Divergence measured the way the projection actually conserves it: from
    // Rhie-Chow face fluxes rather than from wide central differences of
    // cell-centered velocities. See rhieChowFaceU() for why the two differ
    // and which one is the honest number.
    double maxFaceDivergence() const;

private:
    std::size_t index(std::size_t i, std::size_t j) const { return i + j * nx_; }

    // The lid speed at cell-column i: regularized (default) or the
    // classical uniform value, per taperLid_ (see the class comment).
    double lidVelocityAt(std::size_t i) const;

    // Value of `field` at (i+di, j+dj), Dirichlet-mirrored across whichever
    // wall that neighbor would cross (topWallValue is the prescribed value
    // if that wall is the top/lid; every other wall is 0 for both velocity
    // components).
    double dirichletAt(const std::vector<double>& field, std::size_t i, std::size_t j, int di, int dj,
                        double topWallValue) const;

    // Value of `field` at (i+di, j+dj), Neumann-mirrored (ghost equals the
    // boundary cell itself, i.e. zero gradient) across any wall crossed.
    double neumannAt(const std::vector<double>& field, std::size_t i, std::size_t j, int di, int dj) const;

    // **Rhie-Chow interpolation.** The mass flux through the face between
    // cell (i,j) and (i+1,j), which is *not* the plain average of the two
    // cell velocities.
    //
    // Why it is needed, stated precisely for this scheme: the pressure
    // Poisson equation here is solved with the **compact** 5-point
    // Laplacian, (p_E - 2p_P + p_W)/dx^2, but both the divergence that
    // forms its right-hand side and the gradient that corrects the velocity
    // use the **wide** central difference, (p_E - p_W)/(2dx). Composing two
    // wide operators does not reproduce the compact one -- it produces
    // (p_{i+2} - 2p_i + p_{i-2})/(4dx^2), which sees the even and odd
    // sublattices as independent. That mismatch is exactly the classic
    // collocated-grid checkerboard, documented in this class since it was
    // written and finally *measured* by Module 12.2's checkerboardIndex.
    //
    // The Rhie-Chow flux replaces the interpolated pressure-gradient part
    // with the compact face gradient:
    //
    //   u_e = (u_P + u_E)/2 + dt * [ (g_P + g_E)/2 - (p_E - p_P)/dx ]
    //
    // where g is the wide cell-centred gradient the correction step used.
    // Working the divergence of those fluxes through algebraically gives
    // exactly wide_div(u*) - dt * compact_lap(p) -- which is precisely the
    // quantity the Poisson solve drives to zero. So the face divergence of
    // the corrected field is zero to solver tolerance, while the wide
    // cell-centred divergence is not, and never was: it measures a
    // different operator than the one being solved.
    //
    // `dt` is the step size the correction used, kept in lastDt_ because
    // the diagnostic needs it and callers do not pass one. Before the first
    // step it is 0, which makes the correction term vanish and the flux
    // fall back to plain interpolation -- correct, since an unstepped field
    // has no pressure correction to be consistent with.
    // Conservative convective flux of `field` for cell (i,j): the divergence
    // of (u_face * phi_face) over the cell's four faces, with the face value
    // taken from the active scheme. Returns the value that goes on the
    // right-hand side of the predictor, i.e. -div(u phi).
    double conservativeConvection(const std::vector<double>& field, const std::vector<double>& u,
                                   const std::vector<double>& v, std::size_t i, std::size_t j,
                                   double wallValue, double dt) const;
    // One face's value under the active scheme. `upwind`/`downwind` are the
    // two cell values with the flow direction already resolved, `farUpwind`
    // the cell beyond the upwind one -- which a structured grid has and a
    // tetrahedral mesh does not.
    double schemeFaceValue(double upwind, double downwind, double farUpwind) const;

    double rhieChowFaceU(const std::vector<double>& u, std::size_t i, std::size_t j, double dt) const;
    double rhieChowFaceV(const std::vector<double>& v, std::size_t i, std::size_t j, double dt) const;

    // Divergence of cell (i,j) from its four Rhie-Chow face fluxes, with the
    // solid walls contributing exactly zero normal velocity.
    double faceDivergenceAt(const std::vector<double>& u, const std::vector<double>& v, std::size_t i,
                             std::size_t j, double dt) const;

    std::vector<double> applyLaplacian(const std::vector<double>& x) const;
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    void projectToDivergenceFree(std::vector<double>& uStar, std::vector<double>& vStar, double dt);

    std::size_t nx_;
    std::size_t ny_;
    double dx_;
    double dy_;
    double viscosity_;
    double lidVelocity_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> p_;
    double time_ = 0.0;
    ConvectionScheme convection_;
    bool taperLid_;
    double lastDt_ = 0.0; // step size of the most recent correction; see rhieChowFaceU()
};

} // namespace aether::solver
