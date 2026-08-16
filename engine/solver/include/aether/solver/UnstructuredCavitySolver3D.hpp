#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"
#include "aether/solver/UnstructuredFvmBase.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Incompressible Navier-Stokes on an **unstructured tetrahedral mesh** --
// ROADMAP Fase 3, the step where this engine stops being able to simulate
// only boxes.
//
// Fase 2 proved the geometric machinery on the simplest physics available
// (diffusion, against a closed-form Fourier series). That machinery is now
// UnstructuredFvmBase, shared with UnstructuredDiffusionSolver instead of
// copied from it: the pressure Poisson equation *is* the same
// over-relaxed-decomposition Laplacian, built from the same face
// coefficients and the same least-squares gradients. What is new here is
// momentum transport and the projection that couples it to pressure.
//
// **The pressure Poisson equation is now literally that same operator**,
// non-orthogonal deferred correction included. It used to be a second copy
// carrying only the implicit part -- the version Fase 2.2 measured stagnating
// at observed order 0.10 -- which was DIVIDA_TECNICA.md 1.2. Fixing it was a
// change in one place because the base exists.
//
// **Design choices, and why each is what it is:**
//
// - **Collocated storage.** Velocity and pressure both live at cell
//   centroids. A staggered arrangement has no natural analogue on
//   tetrahedra, which is precisely why real unstructured codes are
//   collocated and lean on Rhie-Chow-style face fluxes instead.
//
// - **Upwind convection.** Central differencing of convection is
//   unconditionally unstable above a cell Reynolds number of 2, and this
//   project already measured that biting: the structured cavity at Re=400
//   runs at cell-Re 16.7 with zero margin (see ROADMAP Fase 1). On tetrahedra
//   with irregular cell sizes that margin is worse, not better, so this
//   solver takes the first-order-accurate but bounded choice. Accuracy of
//   the convection scheme is a later concern; a solver that diverges is not
//   a starting point to improve from.
//
// - **Explicit time stepping** for convection, implicit for diffusion, with
//   stableTimeStep() deriving a limit from the actual cell sizes rather than
//   a uniform spacing.
//
// - **Mass conservation is measured on faces, not on a cell-centred
//   difference.** Fase 1 established that the wide cell-centred divergence
//   measures a different operator than the one the projection zeroes; on an
//   unstructured mesh there is no wide stencil to be tempted by, and the
//   face flux is the only meaningful definition. maxFaceDivergence() is the
//   honest diagnostic.
class UnstructuredCavitySolver3D : public UnstructuredFvmBase {
public:
    // Non-orthogonal correctors per projection, rather than iterating the
    // deferred correction to convergence every step: the pressure of an
    // explicit-convection projection scheme is an intermediate quantity, so
    // unstructured codes fix a small number of correctors per step.
    //
    // Measured on the closed cavity, as the peak face divergence the
    // projection leaves behind, with the vortex topology identical in every
    // case:
    //
    //    1 corrector  ->  diverges (2.9e+88)
    //    2            ->  2.284e-02      minimum stable count
    //    4            ->  7.129e-03      the default; +0.8s on a 25s suite
    //    8            ->  8.560e-04
    //   16            ->  1.286e-04
    //
    // The residual falls by roughly half per corrector, which is what says
    // the deferred correction is converging rather than fighting something,
    // so **this diagnostic reports how many correctors were paid for, not a
    // property of the scheme**. One corrector is not merely inaccurate but
    // unstable, which is why the count is not simply minimised.
    //
    // **A caution about reading too much into that table, learned the hard
    // way.** It was first measured on the *channel*, where four correctors
    // left 2.7e-04 of the inflow unaccounted against 7.1e-14 at sixty-four,
    // and that was read as the corrector count being a property of the mesh's
    // non-orthogonality. It was not: the outlet's prescribed pressure was
    // missing from the least-squares gradient stencil, so the correction was
    // being built from a bad gradient exactly where it mattered. With that
    // fixed (DIVIDA_TECNICA.md 2.3) the channel closes to 1e-13 at **two**
    // correctors on every mesh tried, up to non-orthogonality 2.24. What is
    // left above is the closed cavity, which has no Dirichlet face anywhere
    // and pins a reference cell instead -- a genuinely harder problem for
    // this iteration, and the reason the default is four rather than two.
    //
    // Four is therefore the default because it is what the suite is
    // calibrated against; it is **not** a claim that four is enough for a
    // given mesh. lastPressureChange() answers that, and this being a
    // constructor argument is what lets a caller act on the answer.
    static constexpr std::size_t kDefaultPressureCorrectors = 4;

    // `lidVelocity(position)` gives the prescribed wall velocity at a
    // boundary face centroid -- zero on a solid wall, the lid's tangential
    // velocity on the moving one. Passing it as a function keeps the class
    // free of any assumption that the domain is a cube: the same solver runs
    // on an STL-derived mesh by selecting different faces.
    // `isOutlet` marks boundary faces where fluid may *leave* the domain.
    // Without it every boundary is a solid wall: no mass can cross any
    // boundary face, which is fine for a closed cavity and structurally
    // impossible for external flow. An outlet face instead carries a
    // prescribed pressure (Dirichlet, so it enters the Poisson operator's
    // diagonal and right-hand side) and a zero-gradient velocity, letting
    // the projection push mass out through it.
    //
    // An *inlet* needs nothing new: it is a wall with a non-zero prescribed
    // velocity, which `wallVelocity` already expresses.
    // `pressureCorrectors` is the number of non-orthogonal correctors the
    // projection runs per step. It is a constructor argument rather than a
    // constant because **how many are needed is a property of the mesh, not
    // of the solver**, and that was measured rather than assumed -- see
    // kDefaultPressureCorrectors below for the table. lastPressureChange()
    // reports whether the count used was enough.
    //
    // Clamped up to 2: a single corrector does not merely leave a larger
    // residual, it is unstable (the cavity case reaches 2.9e+88), because the
    // correction never catches the pressure it is correcting. Refusing 1
    // outright rather than clamping was the alternative; clamping wins
    // because the only way to ask for 1 is to be guessing, and the failure
    // it produces is a NaN field several thousand steps later rather than an
    // error at the call.
    UnstructuredCavitySolver3D(
        const mesh::TetrahedralMesh& mesh, double viscosity,
        std::function<core::Vector3(const core::Vector3&)> wallVelocity,
        std::function<bool(const core::Vector3&)> isOutlet = {}, double outletPressure = 0.0,
        std::size_t pressureCorrectors = kDefaultPressureCorrectors);

    void step(double dt);

    // **Only the convective limit, because diffusion is implicit.** That is
    // the whole point of the implicit treatment: the explicit viscous bound
    // scales with the square of the smallest cell, and every Delaunay
    // tetrahedralization produces slivers, so it -- not the physics -- was
    // dictating the step. A safety factor is still applied: Fase 1 found the
    // structured cavity at CFL exactly 1.0000 with no margin, which is what
    // made it so fragile.
    double stableTimeStep() const;

    core::Vector3 velocity(std::size_t cell) const { return velocity_.at(cell); }
    double pressure(std::size_t cell) const { return pressure_.at(cell); }
    double time() const { return time_; }

    // Largest |sum of face mass fluxes| / cellVolume. The quantity the
    // projection actually drives to zero -- see the class comment.
    // Includes outlet faces, which do carry flux.
    double maxFaceDivergence() const;

    // Largest per-cell pressure change across the final corrector of the
    // last step. The number that says whether the projection converged or
    // merely ran out of correctors -- the same role lastOuterChange() plays
    // for the diffusion solver, and the reason the corrector count can be
    // chosen from evidence instead of guessed.
    double lastPressureChange() const { return lastPressureChange_; }

    // Net mass flux through every non-wall boundary face: negative where
    // fluid enters, positive where it leaves. For a steady incompressible
    // flow the two must cancel, which is the global check that outlets
    // actually work -- a per-cell divergence can be zero everywhere while
    // the domain as a whole still gains or loses mass.
    double netBoundaryFlux() const;

private:
    // Per-boundary-face condition, indexed alongside the base's
    // boundaryFaces_: the base owns a boundary face's *geometry*, this owns
    // what is prescribed on it.
    struct BoundaryCondition {
        core::Vector3 wallVelocity;
        bool isOutlet = false;
    };

    void buildBoundaryConditions();
    std::vector<double> faceMassFluxes(const std::vector<core::Vector3>& velocity,
                                        const std::vector<double>& pressure, double dt) const;
    void projectToDivergenceFree(std::vector<core::Vector3>& velocityStar, double dt);


    // (V_P/dt + nu * sum_f a_f) x_P - nu * sum_f a_f x_N: the same Laplacian
    // with a shifted diagonal, which is what makes the viscous term implicit.
    // Still symmetric positive definite -- more strongly so than the Poisson
    // operator, since the V/dt shift only adds to the diagonal -- so the same
    // Conjugate Gradient applies unchanged.
    std::vector<double> applyHelmholtzOperator(const std::vector<double>& x, double dt) const;
    std::vector<double> solveHelmholtz(const std::vector<double>& rhs, double dt) const;

    // The boundary mass flux **as the projection left it**, filled in by
    // projectToDivergenceFree(). Stored rather than recomputed because the
    // two are not the same number: the projection corrects an outlet with
    // the compact face gradient a_b (p_outlet - p_P), while recomputing from
    // the cell velocity uses the least-squares gradient instead. Reporting
    // the recomputed one made the domain look like it was losing 13.2% of
    // its inflow -- measuring a different operator than the one being
    // solved, exactly the error ROADMAP Fase 1 diagnosed for the structured
    // cavity's divergence.
    std::vector<double> boundaryFlux_;

    double viscosity_;
    std::function<core::Vector3(const core::Vector3&)> wallVelocity_;
    std::function<bool(const core::Vector3&)> isOutlet_;
    double outletPressure_ = 0.0;

    std::vector<BoundaryCondition> boundaryConditions_;

    std::vector<core::Vector3> velocity_;
    std::vector<double> pressure_;
    double time_ = 0.0;
    std::size_t pressureCorrectors_;
    double lastPressureChange_ = 0.0;
    double lastDt_ = 0.0;
    double maxWallSpeed_ = 0.0;
    bool hasOutlet_ = false;
};

} // namespace aether::solver
