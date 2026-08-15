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
            const double boundaryDistance = d.norm();
            dirichletFaces_.push_back({face.owner, coefficient, boundaryValue_[f],
                                        face.areaVector - d * coefficient, d / boundaryDistance,
                                        boundaryDistance});
            diagonal_[face.owner] += coefficient;
            continue;
        }

        const Vector3 d = mesh_->cellCentroid(face.neighbour) - mesh_->cellCentroid(face.owner);
        const double areaDotD = face.areaVector.dot(d);
        if (areaDotD <= 0.0) {
            continue; // see above: a degenerate/misoriented face contributes nothing
        }
        const double coefficient = face.areaVector.normSquared() / areaDotD;
        // Interpolation weight from the perpendicular distances of the two
        // centroids to the face plane, not a flat 0.5: on a skewed mesh the
        // two cells sit at different distances, and using 0.5 there is a
        // first-order error in every face quantity built from it.
        const double distance = d.norm();
        const Vector3 unitNormal = face.areaVector / face.areaVector.norm();
        const double ownerDistance = std::fabs((face.centroid - mesh_->cellCentroid(face.owner)).dot(unitNormal));
        const double neighbourDistance =
            std::fabs((mesh_->cellCentroid(face.neighbour) - face.centroid).dot(unitNormal));
        const double distanceSum = ownerDistance + neighbourDistance;
        const double ownerWeight = distanceSum > 0.0 ? neighbourDistance / distanceSum : 0.5;
        interiorFaces_.push_back({face.owner, face.neighbour, coefficient,
                                   face.areaVector - d * coefficient, d / distance, distance,
                                   ownerWeight});
        diagonal_[face.owner] += coefficient;
        diagonal_[face.neighbour] += coefficient;
    }

    buildGradientStencils();
}

void UnstructuredDiffusionSolver::buildGradientStencils() {
    const std::size_t n = mesh_->cellCount();
    gradientStencil_.assign(n, {});
    gradientMatrixInverse_.assign(n, {});

    // Accumulate the symmetric normal-equation matrix M = sum_i w_i d_i (x) d_i
    // per cell, alongside the stencil entries the right-hand side needs.
    struct Sym {
        double xx = 0, xy = 0, xz = 0, yy = 0, yz = 0, zz = 0;
    };
    std::vector<Sym> matrices(n);

    const auto accumulate = [](Sym& m, const Vector3& d, double w) {
        m.xx += w * d.x * d.x;
        m.xy += w * d.x * d.y;
        m.xz += w * d.x * d.z;
        m.yy += w * d.y * d.y;
        m.yz += w * d.y * d.z;
        m.zz += w * d.z * d.z;
    };

    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        const auto& face = mesh_->face(f);
        if (mesh_->isBoundaryFace(f)) {
            // Only a Dirichlet face carries information about the gradient;
            // an insulated one asserts the normal derivative is zero, which
            // this least-squares fit would misread as "the value there
            // equals the cell value", biasing the tangential components too.
            if (!boundaryIsDirichlet_[f]) {
                continue;
            }
            const Vector3 d = face.centroid - mesh_->cellCentroid(face.owner);
            const double lengthSquared = d.normSquared();
            if (lengthSquared == 0.0) {
                continue;
            }
            const double w = 1.0 / lengthSquared;
            gradientStencil_[face.owner].push_back({kBoundaryStencil, boundaryValue_[f], d * w});
            accumulate(matrices[face.owner], d, w);
            continue;
        }

        const Vector3 d = mesh_->cellCentroid(face.neighbour) - mesh_->cellCentroid(face.owner);
        const double lengthSquared = d.normSquared();
        if (lengthSquared == 0.0) {
            continue;
        }
        const double w = 1.0 / lengthSquared;
        gradientStencil_[face.owner].push_back({face.neighbour, 0.0, d * w});
        gradientStencil_[face.neighbour].push_back({face.owner, 0.0, -d * w});
        // (-d) (x) (-d) == d (x) d, so both cells accumulate the same term.
        accumulate(matrices[face.owner], d, w);
        accumulate(matrices[face.neighbour], d, w);
    }

    for (std::size_t cell = 0; cell < n; ++cell) {
        const Sym& m = matrices[cell];
        // Cofactors of a symmetric 3x3.
        const double c11 = m.yy * m.zz - m.yz * m.yz;
        const double c12 = m.xz * m.yz - m.xy * m.zz;
        const double c13 = m.xy * m.yz - m.yy * m.xz;
        const double det = m.xx * c11 + m.xy * c12 + m.xz * c13;

        // Rank-deficient when the neighbour directions are too close to
        // coplanar to pin all three components -- scaled by the matrix's own
        // magnitude so the test is dimensionless rather than an absolute
        // threshold that would depend on the mesh's units.
        const double scale = m.xx + m.yy + m.zz;
        if (scale <= 0.0 || std::fabs(det) < 1e-12 * scale * scale * scale) {
            gradientMatrixInverse_[cell].valid = false;
            continue;
        }

        const double inverseDet = 1.0 / det;
        SymmetricInverse& inverse = gradientMatrixInverse_[cell];
        inverse.xx = c11 * inverseDet;
        inverse.xy = c12 * inverseDet;
        inverse.xz = c13 * inverseDet;
        inverse.yy = (m.xx * m.zz - m.xz * m.xz) * inverseDet;
        inverse.yz = (m.xy * m.xz - m.xx * m.yz) * inverseDet;
        inverse.zz = (m.xx * m.yy - m.xy * m.xy) * inverseDet;
        inverse.valid = true;
    }
}

std::vector<Vector3> UnstructuredDiffusionSolver::computeGradients(const std::vector<double>& phi) const {
    std::vector<Vector3> gradients(phi.size());
    std::vector<Vector3> fallback;
    for (std::size_t cell = 0; cell < phi.size(); ++cell) {
        if (!gradientMatrixInverse_[cell].valid) {
            if (fallback.empty()) {
                fallback = computeGradientsGreenGauss(phi);
            }
            gradients[cell] = fallback[cell];
            continue;
        }
        Vector3 rhs;
        for (const GradientStencilEntry& entry : gradientStencil_[cell]) {
            const double neighbourValue =
                entry.neighbour == kBoundaryStencil ? entry.boundaryValue : phi[entry.neighbour];
            rhs += entry.weightedDelta * (neighbourValue - phi[cell]);
        }
        gradients[cell] = gradientMatrixInverse_[cell].apply(rhs);
    }
    return gradients;
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
        const Vector3 averaged =
            gradients[face.owner] * face.ownerWeight + gradients[face.neighbour] * (1.0 - face.ownerWeight);
        const double compactNormalDerivative = (phi[face.neighbour] - phi[face.owner]) / face.distance;
        const Vector3 faceGradient =
            averaged + face.unitD * (compactNormalDerivative - averaged.dot(face.unitD));
        const double flux = faceGradient.dot(face.nonOrthogonalArea);
        correction[face.owner] += flux;
        correction[face.neighbour] -= flux;
    }
    for (const DirichletFace& face : dirichletFaces_) {
        const Vector3& cellGradient = gradients[face.cell];
        const double compactNormalDerivative = (face.value - phi[face.cell]) / face.distance;
        const Vector3 faceGradient =
            cellGradient + face.unitD * (compactNormalDerivative - cellGradient.dot(face.unitD));
        correction[face.cell] += faceGradient.dot(face.nonOrthogonalArea);
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
        lastOuterChange_ = maxChange;
        if (maxChange < tolerance) {
            break;
        }
    }
    return sweep;
}

} // namespace aether::solver
