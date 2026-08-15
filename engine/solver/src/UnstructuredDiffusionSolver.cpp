#include "aether/solver/UnstructuredDiffusionSolver.hpp"

#include <cmath>
#include <stdexcept>

namespace aether::solver {

using core::Vector3;
using mesh::TetrahedralMesh;

UnstructuredDiffusionSolver::UnstructuredDiffusionSolver(const TetrahedralMesh& mesh) : mesh_(&mesh) {
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
    interiorFaces_.clear();
    dirichletFaces_.clear();
    diagonal_.assign(mesh_->cellCount(), 0.0);

    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        const auto& face = mesh_->face(f);

        if (mesh_->isBoundaryFace(f)) {
            if (!boundaryIsDirichlet_[f]) {
                continue; // insulated: no flux, no coefficient, nothing to assemble
            }
            // The "neighbour" of a boundary face is the face centroid
            // itself, where the value is prescribed.
            const Vector3 d = face.centroid - mesh_->cellCentroid(face.owner);
            const double areaDotD = face.areaVector.dot(d);
            if (areaDotD <= 0.0) {
                // Would mean the outward normal points back towards the cell
                // centre, i.e. a broken mesh -- skip rather than emit a
                // negative coefficient that would destroy positive
                // definiteness silently.
                continue;
            }
            const double coefficient = face.areaVector.normSquared() / areaDotD;
            dirichletFaces_.push_back(
                {face.owner, coefficient, boundaryValue_[f], face.areaVector - d * coefficient});
            diagonal_[face.owner] += coefficient;
            continue;
        }

        const Vector3 d = mesh_->cellCentroid(face.neighbour) - mesh_->cellCentroid(face.owner);
        const double areaDotD = face.areaVector.dot(d);
        if (areaDotD <= 0.0) {
            continue; // see above: a degenerate/misoriented face contributes nothing
        }
        const double coefficient = face.areaVector.normSquared() / areaDotD;
        interiorFaces_.push_back(
            {face.owner, face.neighbour, coefficient, face.areaVector - d * coefficient});
        diagonal_[face.owner] += coefficient;
        diagonal_[face.neighbour] += coefficient;
    }
}

std::vector<Vector3> UnstructuredDiffusionSolver::computeGradients(const std::vector<double>& phi) const {
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
        const Vector3 faceGradient = (gradients[face.owner] + gradients[face.neighbour]) * 0.5;
        const double flux = faceGradient.dot(face.nonOrthogonalArea);
        correction[face.owner] += flux;
        correction[face.neighbour] -= flux;
    }
    for (const DirichletFace& face : dirichletFaces_) {
        correction[face.cell] += gradients[face.cell].dot(face.nonOrthogonalArea);
    }
    return correction;
}

double UnstructuredDiffusionSolver::maxNonOrthogonality() const {
    double worst = 0.0;
    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        if (mesh_->isBoundaryFace(f)) {
            continue;
        }
        const auto& face = mesh_->face(f);
        const Vector3 d = mesh_->cellCentroid(face.neighbour) - mesh_->cellCentroid(face.owner);
        const double areaDotD = face.areaVector.dot(d);
        if (areaDotD <= 0.0) {
            continue;
        }
        const Vector3 orthogonal = d * (face.areaVector.normSquared() / areaDotD);
        const Vector3 nonOrthogonal = face.areaVector - orthogonal;
        worst = std::max(worst, nonOrthogonal.norm() / face.areaVector.norm());
    }
    return worst;
}

std::vector<double> UnstructuredDiffusionSolver::applyOperator(const std::vector<double>& x) const {
    // A x, where row P reads
    //   (sum of all face coefficients touching P) * x_P  -  sum_f a_f x_N
    // Symmetric by construction: each interior face contributes the same
    // a_f to both (P,N) and (N,P).
    std::vector<double> result(x.size());
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = diagonal_[cell] * x[cell];
    }
    for (const InteriorFace& face : interiorFaces_) {
        result[face.owner] -= face.coefficient * x[face.neighbour];
        result[face.neighbour] -= face.coefficient * x[face.owner];
    }
    return result;
}

double UnstructuredDiffusionSolver::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

std::size_t UnstructuredDiffusionSolver::solveConjugateGradient(std::size_t maxIterations, double tolerance,
                                                                 std::size_t maxOuterSweeps) {
    if (dirichletFaces_.empty()) {
        throw std::runtime_error(
            "UnstructuredDiffusionSolver: needs at least one Dirichlet boundary face; with pure Neumann "
            "boundaries the operator is singular");
    }

    const std::size_t n = phi_.size();

    std::vector<double> baseRhs(n, 0.0);
    for (const DirichletFace& face : dirichletFaces_) {
        baseRhs[face.cell] += face.coefficient * face.value;
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

        std::vector<double> residual(n);
        {
            const std::vector<double> operatorAppliedToPhi = applyOperator(phi_);
            for (std::size_t i = 0; i < n; ++i) {
                residual[i] = rhs[i] - operatorAppliedToPhi[i];
            }
        }
        std::vector<double> direction = residual;
        double residualDotResidual = dot(residual, residual);

        for (std::size_t iteration = 0; iteration < maxIterations; ++iteration) {
            if (std::sqrt(residualDotResidual) < tolerance) {
                break;
            }
            const std::vector<double> operatorAppliedToDirection = applyOperator(direction);
            const double directionDotOperatorDirection = dot(direction, operatorAppliedToDirection);
            if (directionDotOperatorDirection == 0.0) {
                break;
            }
            const double alpha = residualDotResidual / directionDotOperatorDirection;
            for (std::size_t i = 0; i < n; ++i) {
                phi_[i] += alpha * direction[i];
                residual[i] -= alpha * operatorAppliedToDirection[i];
            }
            const double newResidualDotResidual = dot(residual, residual);
            const double beta = newResidualDotResidual / residualDotResidual;
            for (std::size_t i = 0; i < n; ++i) {
                direction[i] = residual[i] + beta * direction[i];
            }
            residualDotResidual = newResidualDotResidual;
        }

        // Outer convergence: the correction has stopped moving the solution.
        double maxChange = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            maxChange = std::max(maxChange, std::fabs(phi_[i] - previousPhi[i]));
        }
        if (maxChange < tolerance) {
            break;
        }
    }
    return sweep;
}

} // namespace aether::solver
