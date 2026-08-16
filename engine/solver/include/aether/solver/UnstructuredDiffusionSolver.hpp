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
// **The real difficulty: non-orthogonality.** On a Cartesian grid the face
// normal is parallel to the line joining the two cell centres, so the
// diffusive flux through a face is just (phi_N - phi_P)/|d| times the area.
// On tetrahedra it never is, and that misalignment is a property of mesh
// *shape*, not mesh *size* -- refining a tet mesh does not make it more
// orthogonal. The over-relaxed decomposition that handles it now lives in
// UnstructuredFvmBase, which this class shares with the Navier-Stokes
// solver; see that header for the algebra and for why the split was
// extracted.
//
// **The correction is not optional here, and that was measured rather than
// assumed.** A first version assembled only the implicit part. Its mesh
// convergence study stalled: rms error 2.879 -> 2.433 -> 2.365 for n =
// 4/6/8, i.e. observed order 0.42 then 0.10 -- a plateau, not convergence.
// The reason is visible in maxNonOrthogonality() on these meshes (~1.5):
// non-orthogonality is a property of cell *shape*, so refining a tet mesh
// does not reduce it, and the dropped term settles to a mesh-quality error
// floor instead of vanishing. Running that experiment before committing to
// the simpler scheme is the Fase 1 lesson applied.
//
// The correction is handled by **deferred correction**: the implicit part
// keeps the matrix a symmetric positive-definite M-matrix (positive
// diagonal, negative off-diagonals, diagonally dominant once at least one
// Dirichlet face exists) so the base's matrix-free Conjugate Gradient
// applies unchanged, while the non-orthogonal term is evaluated from the
// previous iterate and moved to the right-hand side. Outer sweeps repeat
// until it stops changing. The face gradient it needs comes from
// least-squares cell gradients averaged to the face.
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
    std::vector<core::Vector3> cellGradients() const;

    // Largest per-cell change across the final outer sweep. The number that
    // says whether the deferred-correction loop actually converged or merely
    // ran out of sweeps -- without it, a solve that stopped at the iteration
    // cap is indistinguishable from one that settled, and the two mean very
    // different things about the result.
    double lastOuterChange() const { return lastOuterChange_; }

    double value(std::size_t cell) const { return phi_.at(cell); }

private:
    std::vector<double> applyOperator(const std::vector<double>& x) const;
    void rebuildCoefficients();

    bool isDirichlet(const BoundaryFace& face) const { return boundaryIsDirichlet_[face.meshFace]; }
    double dirichletValue(const BoundaryFace& face) const { return boundaryValue_[face.meshFace]; }

    // The non-orthogonal flux the implicit part cannot represent, evaluated
    // from the previous iterate and accumulated per cell.
    std::vector<double> nonOrthogonalCorrection(const std::vector<double>& phi) const;

    // Green-Gauss: cheap, but only first-order accurate on a skewed mesh --
    // kept as the fallback for the rare cell whose least-squares stencil is
    // rank-deficient (too few, or nearly coplanar, neighbours). Returning a
    // zero gradient there instead would be a wrong answer that looks like a
    // result; see UnstructuredFvmBase::leastSquaresGradients.
    std::vector<core::Vector3> computeGradientsGreenGauss(const std::vector<double>& phi) const;

    // Least-squares cell gradients, with the Green-Gauss fallback supplied
    // only when this mesh actually has a deficient stencil to fall back for.
    std::vector<core::Vector3> computeGradients(const std::vector<double>& phi) const;

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
