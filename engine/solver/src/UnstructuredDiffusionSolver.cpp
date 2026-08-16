#include "aether/solver/UnstructuredDiffusionSolver.hpp"

#include <cmath>
#include <stdexcept>

namespace aether::solver {

using core::Vector3;
using mesh::TetrahedralMesh;

UnstructuredDiffusionSolver::UnstructuredDiffusionSolver(const TetrahedralMesh& mesh)
    : UnstructuredFvmBase(mesh) {
    phi_.assign(mesh.cellCount(), 0.0);
    boundaryIsDirichlet_.assign(mesh.faceCount(), false);
    boundaryValue_.assign(mesh.faceCount(), 0.0);
    rebuildCoefficients();
}

void UnstructuredDiffusionSolver::setDirichletBoundary(
    const std::function<bool(const Vector3&)>& selector, double value) {
    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        if (!mesh_->isBoundaryFace(f)) {
            continue;
        }
        if (selector(mesh_->face(f).centroid)) {
            boundaryIsDirichlet_[f] = true;
            boundaryValue_[f] = value;
        }
    }
    rebuildCoefficients();
}

void UnstructuredDiffusionSolver::rebuildCoefficients() {
    // A Dirichlet face is an extra connection and enters the operator; an
    // insulated one carries no flux, so it contributes no coefficient at all.
    buildFaceGeometry([this](std::size_t meshFace) { return boundaryIsDirichlet_[meshFace]; });

    hasDirichletFace_ = false;
    for (const BoundaryFace& face : boundaryFaces_) {
        if (isDirichlet(face)) {
            hasDirichletFace_ = true;
            break;
        }
    }

    buildGradientStencils([this](std::size_t meshFace, double& value) {
        if (!boundaryIsDirichlet_[meshFace]) {
            return false;
        }
        value = boundaryValue_[meshFace];
        return true;
    });
}

std::vector<Vector3> UnstructuredDiffusionSolver::computeGradients(const std::vector<double>& phi) const {
    if (!hasDeficientStencil()) {
        return leastSquaresGradients(phi);
    }
    const std::vector<Vector3> fallback = computeGradientsGreenGauss(phi);
    return leastSquaresGradients(phi, &fallback);
}

std::vector<Vector3> UnstructuredDiffusionSolver::computeGradientsGreenGauss(
    const std::vector<double>& phi) const {
    // Green-Gauss: grad(phi)_P = (1/V_P) * sum_f phi_f A_f_out. The face
    // value is a plain average for interior faces, the prescribed value on a
    // Dirichlet face, and phi_P itself on an insulated face (zero gradient
    // means the face carries the cell's own value).
    std::vector<Vector3> gradients(phi.size());
    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        const auto& face = mesh_->face(f);
        if (mesh_->isBoundaryFace(f)) {
            const double faceValue = boundaryIsDirichlet_[f] ? boundaryValue_[f] : phi[face.owner];
            gradients[face.owner] += face.areaVector * faceValue;
            continue;
        }
        const double faceValue = 0.5 * (phi[face.owner] + phi[face.neighbour]);
        gradients[face.owner] += face.areaVector * faceValue;
        gradients[face.neighbour] -= face.areaVector * faceValue;
    }
    for (std::size_t cell = 0; cell < gradients.size(); ++cell) {
        gradients[cell] = gradients[cell] / mesh_->cellVolume(cell);
    }
    return gradients;
}

std::vector<Vector3> UnstructuredDiffusionSolver::cellGradients() const { return computeGradients(phi_); }

std::vector<double> UnstructuredDiffusionSolver::nonOrthogonalCorrection(
    const std::vector<double>& phi) const {
    // Evaluated from the previous iterate and moved to the right-hand side,
    // which is what keeps the matrix itself symmetric and positive definite.
    const std::vector<Vector3> gradients = computeGradients(phi);

    std::vector<double> correction(phi.size(), 0.0);
    for (const InteriorFace& face : interiorFaces_) {
        // **Corrected face gradient.** The plain average approximates the
        // gradient at the midpoint of the centroid segment, not at the face,
        // which is an O(h) error on a skewed mesh -- and since
        // A_nonorth is not perpendicular to d under the over-relaxed
        // decomposition, that error leaks straight into the flux. The fix is
        // to keep the averaged *tangential* part but replace the component
        // along d by the compact difference (phi_N - phi_P)/|d|, which is the
        // one direction the two cell values pin down accurately. Standard
        // practice, and the remaining first-order term identified when the
        // least-squares gradient alone left the order at ~1.
        const Vector3 averaged = averagedFaceGradient(face, gradients);
        const double compactNormalDerivative = (phi[face.neighbour] - phi[face.owner]) / face.distance;
        const Vector3 faceGradient =
            averaged + face.unitD * (compactNormalDerivative - averaged.dot(face.unitD));
        const double flux = faceGradient.dot(face.nonOrthogonalArea);
        correction[face.owner] += flux;
        correction[face.neighbour] -= flux;
    }
    for (const BoundaryFace& face : boundaryFaces_) {
        if (!isDirichlet(face)) {
            continue;
        }
        const Vector3& cellGradient = gradients[face.cell];
        const double compactNormalDerivative = (dirichletValue(face) - phi[face.cell]) / face.distance;
        const Vector3 faceGradient =
            cellGradient + face.unitD * (compactNormalDerivative - cellGradient.dot(face.unitD));
        correction[face.cell] += faceGradient.dot(face.nonOrthogonalArea);
    }
    return correction;
}

std::vector<double> UnstructuredDiffusionSolver::applyOperator(const std::vector<double>& x) const {
    // A x, where row P reads
    //   (sum of all face coefficients touching P) * x_P  -  sum_f a_f x_N
    // Symmetric by construction: each interior face contributes the same
    // a_f to both (P,N) and (N,P).
    std::vector<double> result(x.size());
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = laplacianDiagonal_[cell] * x[cell];
    }
    for (const InteriorFace& face : interiorFaces_) {
        result[face.owner] -= face.coefficient * x[face.neighbour];
        result[face.neighbour] -= face.coefficient * x[face.owner];
    }
    return result;
}

std::size_t UnstructuredDiffusionSolver::solveConjugateGradient(std::size_t maxIterations, double tolerance,
                                                                 std::size_t maxOuterSweeps) {
    if (!hasDirichletFace_) {
        throw std::runtime_error(
            "UnstructuredDiffusionSolver: needs at least one Dirichlet boundary face; with pure Neumann "
            "boundaries the operator is singular");
    }

    const std::size_t n = phi_.size();

    std::vector<double> baseRhs(n, 0.0);
    for (const BoundaryFace& face : boundaryFaces_) {
        if (isDirichlet(face)) {
            baseRhs[face.cell] += face.coefficient * dirichletValue(face);
        }
    }

    std::size_t sweep = 0;
    for (; sweep < maxOuterSweeps; ++sweep) {
        // Deferred correction: the non-orthogonal flux is evaluated from the
        // current iterate and treated as a source, so the matrix stays SPD
        // and plain CG still applies. On the first sweep phi_ is whatever it
        // was (zero for a fresh solver), which makes that sweep the
        // orthogonal-only solution and every later one a refinement of it.
        const std::vector<double> correction = nonOrthogonalCorrection(phi_);
        std::vector<double> rhs(n);
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] = baseRhs[i] + correction[i];
        }

        const std::vector<double> previousPhi = phi_;

        conjugateGradient([this](const std::vector<double>& v) { return applyOperator(v); }, rhs, phi_,
                          maxIterations, tolerance);

        // Outer convergence: the correction has stopped moving the solution.
        double maxChange = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            maxChange = std::max(maxChange, std::fabs(phi_[i] - previousPhi[i]));
        }
        lastOuterChange_ = maxChange;
        if (maxChange < tolerance) {
            break;
        }
    }
    return sweep;
}

} // namespace aether::solver
