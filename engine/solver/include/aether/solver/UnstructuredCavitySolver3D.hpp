#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Incompressible Navier-Stokes on an **unstructured tetrahedral mesh** --
// ROADMAP Fase 3, the step where this engine stops being able to simulate
// only boxes.
//
// Fase 2 proved the geometric machinery on the simplest physics available
// (diffusion, against a closed-form Fourier series). This class reuses that
// machinery unchanged for the hardest part of an incompressible solve: the
// pressure Poisson equation *is* the same over-relaxed-decomposition
// Laplacian, with the same non-orthogonal deferred correction and the same
// least-squares gradients. What is new here is momentum transport and the
// projection that couples it to pressure.
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
// - **Explicit time stepping**, matching every other Navier-Stokes solver
//   here, with stableTimeStep() deriving a limit from the actual cell sizes
//   rather than a uniform spacing.
//
// - **Mass conservation is measured on faces, not on a cell-centred
//   difference.** Fase 1 established that the wide cell-centred divergence
//   measures a different operator than the one the projection zeroes; on an
//   unstructured mesh there is no wide stencil to be tempted by, and the
//   face flux is the only meaningful definition. maxFaceDivergence() is the
//   honest diagnostic.
class UnstructuredCavitySolver3D {
public:
    // `lidVelocity(position)` gives the prescribed wall velocity at a
    // boundary face centroid -- zero on a solid wall, the lid's tangential
    // velocity on the moving one. Passing it as a function keeps the class
    // free of any assumption that the domain is a cube: the same solver runs
    // on an STL-derived mesh by selecting different faces.
    UnstructuredCavitySolver3D(const mesh::TetrahedralMesh& mesh, double viscosity,
                                std::function<core::Vector3(const core::Vector3&)> wallVelocity);

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
    std::size_t cellCount() const { return velocity_.size(); }

    // Largest |sum of face mass fluxes| / cellVolume. The quantity the
    // projection actually drives to zero -- see the class comment.
    double maxFaceDivergence() const;

private:
    struct InteriorFace {
        std::size_t owner;
        std::size_t neighbour;
        double laplacianCoefficient; // a_f = |A|^2 / (A . d)
        core::Vector3 areaVector;    // owner -> neighbour
        core::Vector3 nonOrthogonalArea;
        core::Vector3 unitD;
        double distance;
        double ownerWeight;
    };

    struct BoundaryFace {
        std::size_t cell;
        core::Vector3 areaVector; // outward
        core::Vector3 centroid;
        core::Vector3 wallVelocity;
        double laplacianCoefficient; // for the viscous wall flux
        double distance;
    };

    void buildFaces();
    void buildGradientStencils();
    std::vector<core::Vector3> scalarGradients(const std::vector<double>& field) const;
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

    const mesh::TetrahedralMesh* mesh_;
    double viscosity_;
    std::function<core::Vector3(const core::Vector3&)> wallVelocity_;

    std::vector<InteriorFace> interiorFaces_;
    std::vector<BoundaryFace> boundaryFaces_;
    std::vector<double> poissonDiagonal_;

    struct GradientStencilEntry {
        std::size_t neighbour;
        core::Vector3 weightedDelta;
    };
    struct SymmetricInverse {
        double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
        bool valid = false;
        core::Vector3 apply(const core::Vector3& v) const {
            return {xx * v.x + xy * v.y + xz * v.z, xy * v.x + yy * v.y + yz * v.z,
                    xz * v.x + yz * v.y + zz * v.z};
        }
    };
    std::vector<std::vector<GradientStencilEntry>> gradientStencil_;
    std::vector<SymmetricInverse> gradientMatrixInverse_;

    std::vector<core::Vector3> velocity_;
    std::vector<double> pressure_;
    double time_ = 0.0;
    double lastDt_ = 0.0;
    double maxWallSpeed_ = 0.0;
};

} // namespace aether::solver
