#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Shared finite-volume machinery for every solver in this project that runs
// on a tetrahedral mesh, extracted from UnstructuredDiffusionSolver and
// UnstructuredCavitySolver3D.
//
// **Why it was extracted, and why now.** The two classes had grown
// independent copies of the same five things: the over-relaxed face
// decomposition, the face coefficient |A|^2/(A.d), the inverse-distance
// least-squares gradient stencil, the symmetric 3x3 inverse with its rank
// guard, and the matrix-free Conjugate Gradient loop. That is the same
// situation that forced StaggeredCavityBase3D out of the six 3D cavities,
// and the reason recorded there applies unchanged: **porting a correction N
// times is the expensive version of the problem.** It was not a
// hypothetical risk here -- the copies had already diverged within a day of
// being written, in two places, both recorded in DIVIDA_TECNICA.md:
//
//   1.2  the diffusion solver's Laplacian carries the non-orthogonal
//        deferred correction; the Navier-Stokes pressure Poisson operator,
//        which is the *same operator*, does not -- and Fase 2.2 measured
//        that omission making the error stagnate at order 0.10.
//   1.3  a rank-deficient gradient stencil falls back to Green-Gauss in one
//        copy and silently returns a zero gradient in the other.
//
// Neither is fixed here: this class is what makes fixing them *once*
// possible, which is why the roadmap ordered it first. Where the two copies
// genuinely differ, the difference is now a named argument at the call site
// (see buildFaceGeometry and buildGradientStencils) rather than a silent
// divergence between two bodies of code nobody diffs.
//
// **The extraction is behavior-preserving, and that was designed for rather
// than hoped for.** Every shared loop below visits mesh faces in the same
// order the two originals did and accumulates into the same sums in the same
// sequence, so the floating-point results are bit-for-bit what they were.
// That is why buildFaceGeometry takes a predicate instead of leaving each
// solver to add its own boundary contributions afterwards: assembling
// interiors first and boundaries second would be equally correct and would
// have perturbed every sum in the last bit, which destroys the only cheap
// way to check that a refactor changed nothing.
//
// **No virtual functions**, for the same reason StaggeredCavityBase3D has
// none: the customisation points are build-time (std::function, called once
// per face while assembling) and the hot path is a template
// (conjugateGradient), so nothing in an inner loop pays for dispatch. The
// destructor is protected and non-virtual -- this is a base to reuse code,
// not to delete through.
class UnstructuredFvmBase {
public:
    std::size_t cellCount() const { return mesh_->cellCount(); }

    // Largest |A_nonorth| / |A| over the interior faces: 0 on a perfectly
    // orthogonal mesh, approaching 1 as a face becomes parallel to the line
    // between the cell centres. Exposed because it is the quantity that
    // decides whether the non-orthogonal correction matters, and a caller
    // measuring accuracy should be able to report it alongside.
    double maxNonOrthogonality() const;

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
    // operator, which is the predicate below.
    struct BoundaryFace {
        std::size_t meshFace; // index into the mesh's face list
        std::size_t cell;
        double coefficient;              // a_b = |A|^2 / (A . d_b)
        core::Vector3 areaVector;        // outward
        core::Vector3 centroid;
        core::Vector3 nonOrthogonalArea;
        core::Vector3 unitD;
        double distance;
    };

    // Fills interiorFaces_, boundaryFaces_ and laplacianDiagonal_ in a single
    // pass over the mesh.
    //
    // `boundaryEntersOperator(meshFace)` decides whether a boundary face
    // contributes its coefficient to the diagonal. That is the one place the
    // two solvers' assembly genuinely differs, and both readings are the same
    // rule stated for a different condition: a **prescribed value**
    // (Dirichlet -- a fixed temperature, an outlet pressure) makes the face an
    // extra connection and belongs in the operator, while a **zero-gradient**
    // face (an insulated wall, the pressure condition on a solid wall)
    // carries no flux and therefore no coefficient at all. Passing {} means
    // no boundary face enters, which leaves the operator singular.
    //
    // Faces with A . d <= 0 are skipped entirely rather than assembled with a
    // negative coefficient: that would mean the outward normal points back
    // towards the cell centre, i.e. a broken mesh, and a negative coefficient
    // destroys positive definiteness silently -- the CG below would then
    // converge to something plausible and wrong.
    void buildFaceGeometry(const std::function<bool(std::size_t)>& boundaryEntersOperator);

    // -- Least-squares gradients ------------------------------------------

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

    // Builds the inverse-distance-weighted least-squares stencil: the fit
    // that best reproduces the measured differences to every neighbour,
    // minimising sum_i w_i (grad(phi)_P . d_i - (phi_i - phi_P))^2 with
    // w_i = 1/|d_i|^2. Unlike Green-Gauss it stays second-order on a skewed
    // mesh, which is why it exists: the Green-Gauss version's convergence
    // order stalled around 0.6-1.0 instead of reaching 2 (Fase 2.2).
    //
    // `prescribedBoundaryValue(meshFace, value)` returns true, writing
    // `value`, for a boundary face whose value is known. Only such a face
    // carries information about the gradient -- a zero-gradient face asserts
    // that the *normal derivative* vanishes, which this fit would misread as
    // "the value there equals the cell value", biasing the tangential
    // components too. Pass {} for an interior-only stencil.
    void buildGradientStencils(
        const std::function<bool(std::size_t, double&)>& prescribedBoundaryValue = {});

    // Gradient per cell from the least-squares fit. Cells whose stencil is
    // rank-deficient take `fallback` when one is supplied, and are left at
    // zero when it is not.
    //
    // **Leaving them at zero is not a neutral choice**, and a caller that
    // does so is making it explicitly for that reason: a zero gradient is not
    // "no answer", it is a wrong answer that looks plausible, and deficient
    // stencils appear exactly where the geometry is complicated and the
    // gradient matters most (DIVIDA_TECNICA.md 1.3). hasDeficientStencil()
    // says whether the question arises at all on this mesh, so the fallback
    // is only computed when it does.
    std::vector<core::Vector3> leastSquaresGradients(
        const std::vector<double>& field, const std::vector<core::Vector3>* fallback = nullptr) const;

    bool hasDeficientStencil() const { return hasDeficientStencil_; }

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

    // -- Linear solve ------------------------------------------------------

    static double dot(const std::vector<double>& a, const std::vector<double>& b);

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

    const mesh::TetrahedralMesh* mesh_;
    std::vector<InteriorFace> interiorFaces_;
    std::vector<BoundaryFace> boundaryFaces_;
    // sum_f a_f per cell: the diagonal of the assembled Laplacian, over the
    // interior faces plus whichever boundary faces enter the operator.
    std::vector<double> laplacianDiagonal_;
    std::vector<std::vector<GradientStencilEntry>> gradientStencil_;
    std::vector<SymmetricInverse> gradientMatrixInverse_;
    bool hasDeficientStencil_ = false;
};

} // namespace aether::solver
