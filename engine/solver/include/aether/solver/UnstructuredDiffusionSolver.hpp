#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"

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
// orthogonal. This class uses the standard over-relaxed decomposition,
// splitting the face area vector A into a part along d = c_N - c_P and a
// remainder:
//
//   A_orth = (|A|^2 / (A . d)) d,      A_nonorth = A - A_orth
//
// so that the flux is
//
//   (grad phi)_f . A  =  a_f (phi_N - phi_P)  +  (grad phi)_f . A_nonorth
//                        \___ implicit ____/     \___ correction ____/
//
// with a_f = |A|^2 / (A . d) > 0.
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
// Dirichlet face exists) so the project's existing matrix-free Conjugate
// Gradient applies unchanged, while the non-orthogonal term is evaluated
// from the previous iterate and moved to the right-hand side. Outer sweeps
// repeat until it stops changing. The face gradient it needs comes from
// Green-Gauss cell gradients averaged to the face.
class UnstructuredDiffusionSolver {
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
    std::size_t cellCount() const { return phi_.size(); }

    // Largest |A_nonorth| / |A| over the interior faces: 0 on a perfectly
    // orthogonal mesh, approaching 1 as a face becomes parallel to the line
    // between the cell centres. Exposed because it is the quantity that
    // decides whether the omitted correction matters, and a caller
    // measuring accuracy should be able to report it alongside.
    double maxNonOrthogonality() const;

private:
    std::vector<double> applyOperator(const std::vector<double>& x) const;
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    void rebuildCoefficients();
    void buildGradientStencils();

    // The non-orthogonal flux the implicit part cannot represent, evaluated
    // from the previous iterate and accumulated per cell.
    std::vector<double> nonOrthogonalCorrection(const std::vector<double>& phi) const;

    // Green-Gauss: cheap, but only first-order accurate on a skewed mesh --
    // kept as the fallback for the rare cell whose least-squares stencil is
    // rank-deficient (too few, or nearly coplanar, neighbours).
    std::vector<core::Vector3> computeGradientsGreenGauss(const std::vector<double>& phi) const;

    // Inverse-distance-weighted least squares: fits the gradient that best
    // reproduces the measured differences to every neighbour, minimising
    // sum_i w_i (grad(phi)_P . d_i - (phi_i - phi_P))^2 with w_i = 1/|d_i|^2.
    // Unlike Green-Gauss it stays second-order on a skewed mesh, which is
    // why it exists here: the Green-Gauss version's convergence order
    // stalled around 0.6-1.0 instead of reaching 2.
    std::vector<core::Vector3> computeGradients(const std::vector<double>& phi) const;

    struct GradientStencilEntry {
        std::size_t neighbour; // kBoundaryStencil when the "neighbour" is a Dirichlet face
        double boundaryValue;
        core::Vector3 weightedDelta; // w_i * d_i
    };
    static constexpr std::size_t kBoundaryStencil = static_cast<std::size_t>(-1);

    // Inverse of the symmetric 3x3 normal-equation matrix, stored as its six
    // unique components (xx, xy, xz, yy, yz, zz). Geometry only, so it is
    // built once alongside the face coefficients.
    struct SymmetricInverse {
        double xx = 0.0, xy = 0.0, xz = 0.0, yy = 0.0, yz = 0.0, zz = 0.0;
        bool valid = false; // false => stencil was rank-deficient, use Green-Gauss for this cell
        core::Vector3 apply(const core::Vector3& v) const {
            return {xx * v.x + xy * v.y + xz * v.z, xy * v.x + yy * v.y + yz * v.z,
                    xz * v.x + yz * v.y + zz * v.z};
        }
    };

    struct InteriorFace {
        std::size_t owner;
        std::size_t neighbour;
        double coefficient;              // a_f = |A|^2 / (A . d)
        core::Vector3 nonOrthogonalArea; // A - (|A|^2/(A.d)) d, as seen from the owner
        core::Vector3 unitD;             // normalised c_N - c_P
        double distance;                 // |c_N - c_P|
        double ownerWeight;              // distance-based interpolation weight for the owner
    };

    struct DirichletFace {
        std::size_t cell;
        double coefficient; // |A|^2 / (A . d_b), with d_b = faceCentroid - cellCentroid
        double value;
        core::Vector3 nonOrthogonalArea;
        core::Vector3 unitD;
        double distance;
    };

    const mesh::TetrahedralMesh* mesh_;
    std::vector<InteriorFace> interiorFaces_;
    std::vector<DirichletFace> dirichletFaces_;
    // Per-boundary-face prescribed value, or "no value" when insulated.
    std::vector<bool> boundaryIsDirichlet_;
    std::vector<double> boundaryValue_;
    std::vector<double> diagonal_;
    std::vector<double> phi_;
    double lastOuterChange_ = 0.0;
    std::vector<std::vector<GradientStencilEntry>> gradientStencil_;
    std::vector<SymmetricInverse> gradientMatrixInverse_;
};

} // namespace aether::solver
