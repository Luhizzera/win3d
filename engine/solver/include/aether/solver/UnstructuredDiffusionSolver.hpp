#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"
#include "aether/solver/UnstructuredFvmBase.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Steady diffusion (Laplace's equation) by the finite-volume method on an
// **unstructured tetrahedral mesh** -- the first solver in this project to
// consume mesh generation rather than a uniform Cartesian grid.
//
// **Why this class is the point of ROADMAP Fase 2.** Every other solver
// here takes StructuredGrid3D, so until now the engine could import an STL,
// tetrahedralize it, and then solve in a box. This closes that gap for the
// simplest physics the project has, deliberately: diffusion already has an
// independent analytical answer (the 2D Fourier series that
// SteadyDiffusionSolver is validated against), so solving the *same*
// problem on tetrahedra isolates what changed -- the geometric
// discretization -- from what did not -- the physics.
//
// **What is left in this class is only the physics.** The geometry (the
// over-relaxed decomposition, the non-orthogonal deferred correction, the
// least-squares gradients with their Green-Gauss fallback, the Laplacian and
// the Conjugate Gradient that solves it) is UnstructuredFvmBase, shared with
// the Navier-Stokes solver. What remains here is the boundary condition --
// which faces are held at a fixed value and which are insulated -- and the
// right-hand side that follows from it. That split is the whole reason the
// base exists; see its header.
//
// **The non-orthogonality this rests on.** On a Cartesian grid the face
// normal is parallel to the line joining the two cell centres, so the
// diffusive flux through a face is just (phi_N - phi_P)/|d| times the area.
// On tetrahedra it never is, and that misalignment is a property of mesh
// *shape*, not mesh *size* -- refining a tet mesh does not make it more
// orthogonal (maxNonOrthogonality() stays ~1.5 on these meshes). Fase 2.2
// measured what dropping the resulting correction costs: rms error 2.879 ->
// 2.433 -> 2.365 for n = 4/6/8, an observed order of 0.42 then 0.10 -- a
// plateau, not convergence. Running that experiment before committing to the
// simpler scheme is the Fase 1 lesson applied.
class UnstructuredDiffusionSolver : public UnstructuredFvmBase {
public:
    explicit UnstructuredDiffusionSolver(const mesh::TetrahedralMesh& mesh);

    // Every boundary face starts insulated (zero-gradient / Neumann), which
    // contributes no flux and therefore needs no coefficient at all. Faces
    // whose *centroid* satisfies `selector` are instead fixed to `value`.
    // Selecting by centroid rather than by index is what lets a caller say
    // "the face on the x = 0 plane" without knowing anything about how the
    // mesh generator happened to number things. A later call overrides an
    // earlier one for the same face.
    void setDirichletBoundary(const std::function<bool(const core::Vector3&)>& selector, double value);

    // Outer deferred-correction sweeps, each an inner matrix-free CG solve
    // (same pattern as SteadyDiffusionSolver). Returns the number of outer
    // sweeps actually run; sweeps stop early once the solution stops moving
    // by more than `tolerance` between them.
    //
    // Requires at least one Dirichlet face: with pure Neumann boundaries the
    // operator is singular (the solution is defined only up to a constant),
    // and this class does not pin a reference cell the way the pressure
    // solvers do -- the diffusion problems it is meant for always have a
    // fixed value somewhere.
    std::size_t solveConjugateGradient(std::size_t maxIterations = 20000, double tolerance = 1e-10,
                                        std::size_t maxOuterSweeps = 50);

    // Gradient at each cell from the current solution. Exposed because it is
    // the quantity the non-orthogonal correction is built from, so a caller
    // checking this class's accuracy can inspect it.
    std::vector<core::Vector3> cellGradients() const { return computeCellGradients(phi_); }

    // Largest per-cell change across the final outer sweep. The number that
    // says whether the deferred-correction loop actually converged or merely
    // ran out of sweeps -- without it, a solve that stopped at the iteration
    // cap is indistinguishable from one that settled, and the two mean very
    // different things about the result.
    double lastOuterChange() const { return lastOuterChange_; }

    double value(std::size_t cell) const { return phi_.at(cell); }

private:
    void rebuildCoefficients();

    // Per-boundary-face prescribed value, or "no value" when insulated.
    // Indexed by *mesh* face, since that is what a selector resolves to and
    // what the base hands back in BoundaryFace::meshFace.
    std::vector<bool> boundaryIsDirichlet_;
    std::vector<double> boundaryValue_;
    bool hasDirichletFace_ = false;
    std::vector<double> phi_;
    double lastOuterChange_ = 0.0;
};

} // namespace aether::solver
