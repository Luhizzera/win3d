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
// **One divergence from the shared base is still real and is not fixed
// here**: applyPoissonOperator() assembles only the implicit part of that
// Laplacian, without the non-orthogonal deferred correction the diffusion
// solver applies to the same operator -- the correction Fase 2.2 measured as
// the difference between converging and stagnating at order 0.10. It is
// recorded as DIVIDA_TECNICA.md 1.2, and sharing the base is what makes
// fixing it a change in one place rather than two.
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
    UnstructuredCavitySolver3D(
        const mesh::TetrahedralMesh& mesh, double viscosity,
        std::function<core::Vector3(const core::Vector3&)> wallVelocity,
        std::function<bool(const core::Vector3&)> isOutlet = {}, double outletPressure = 0.0);

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
    std::vector<double> applyPoissonOperator(const std::vector<double>& x) const;

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
    double lastDt_ = 0.0;
    double maxWallSpeed_ = 0.0;
    bool hasOutlet_ = false;
};

} // namespace aether::solver
