#include "aether/solver/UnstructuredFvmBase.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

using core::Vector3;

void UnstructuredFvmBase::buildFaceGeometry(
    const std::function<bool(std::size_t)>& boundaryEntersOperator) {
    interiorFaces_.clear();
    boundaryFaces_.clear();
    interiorDiagonal_.assign(mesh_->cellCount(), 0.0);
    laplacianDiagonal_.assign(mesh_->cellCount(), 0.0);

    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        const auto& face = mesh_->face(f);

        if (mesh_->isBoundaryFace(f)) {
            // d_b runs from the cell centroid to the face centroid, where a
            // prescribed value would live.
            const Vector3 d = face.centroid - mesh_->cellCentroid(face.owner);
            const double areaDotD = face.areaVector.dot(d);
            if (areaDotD <= 0.0) {
                continue; // misoriented or degenerate: see the header
            }
            const double coefficient = face.areaVector.normSquared() / areaDotD;
            const double distance = d.norm();
            BoundaryFace boundary;
            boundary.meshFace = f;
            boundary.cell = face.owner;
            boundary.coefficient = coefficient;
            boundary.areaVector = face.areaVector;
            boundary.centroid = face.centroid;
            boundary.nonOrthogonalArea = face.areaVector - d * coefficient;
            boundary.unitD = d / distance;
            boundary.distance = distance;
            boundary.entersOperator = boundaryEntersOperator && boundaryEntersOperator(f);
            boundaryFaces_.push_back(boundary);
            if (boundary.entersOperator) {
                laplacianDiagonal_[face.owner] += coefficient;
            }
            continue;
        }

        const Vector3 d = mesh_->cellCentroid(face.neighbour) - mesh_->cellCentroid(face.owner);
        const double areaDotD = face.areaVector.dot(d);
        if (areaDotD <= 0.0) {
            continue;
        }
        const double coefficient = face.areaVector.normSquared() / areaDotD;
        const double distance = d.norm();

        InteriorFace interior;
        interior.owner = face.owner;
        interior.neighbour = face.neighbour;
        interior.coefficient = coefficient;
        interior.areaVector = face.areaVector;
        interior.nonOrthogonalArea = face.areaVector - d * coefficient;
        interior.unitD = d / distance;
        interior.distance = distance;

        const Vector3 unitNormal = face.areaVector / face.areaVector.norm();
        const double ownerDistance =
            std::fabs((face.centroid - mesh_->cellCentroid(face.owner)).dot(unitNormal));
        const double neighbourDistance =
            std::fabs((mesh_->cellCentroid(face.neighbour) - face.centroid).dot(unitNormal));
        const double distanceSum = ownerDistance + neighbourDistance;
        interior.ownerWeight = distanceSum > 0.0 ? neighbourDistance / distanceSum : 0.5;

        interiorDiagonal_[face.owner] += coefficient;
        interiorDiagonal_[face.neighbour] += coefficient;
        laplacianDiagonal_[face.owner] += coefficient;
        laplacianDiagonal_[face.neighbour] += coefficient;
        interiorFaces_.push_back(interior);
    }
}

void UnstructuredFvmBase::buildGradientStencils(bool useBoundaryValues) {
    const std::size_t n = mesh_->cellCount();
    gradientStencil_.assign(n, {});
    gradientMatrixInverse_.assign(n, {});
    deficientStencilCount_ = 0;

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
            double prescribed = 0.0;
            if (!useBoundaryValues || !boundaryFaceValue_ || !boundaryFaceValue_(f, prescribed)) {
                continue; // zero-gradient face: carries no gradient information
            }
            const Vector3 d = face.centroid - mesh_->cellCentroid(face.owner);
            const double lengthSquared = d.normSquared();
            if (lengthSquared == 0.0) {
                continue;
            }
            const double w = 1.0 / lengthSquared;
            gradientStencil_[face.owner].push_back({kBoundaryStencil, prescribed, d * w});
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
            ++deficientStencilCount_;
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

std::vector<Vector3> UnstructuredFvmBase::computeCellGradients(const std::vector<double>& field) const {
    // The Green-Gauss fallback is only built when this mesh actually has a
    // cell that needs it -- it costs a full pass over the faces.
    std::vector<Vector3> fallback;
    if (deficientStencilCount_ > 0) {
        fallback = greenGaussGradients(field);
    }

    std::vector<Vector3> gradients(field.size());
    for (std::size_t cell = 0; cell < field.size(); ++cell) {
        if (!gradientMatrixInverse_[cell].valid) {
            gradients[cell] = fallback[cell];
            continue;
        }
        Vector3 rhs;
        for (const GradientStencilEntry& entry : gradientStencil_[cell]) {
            const double neighbourValue =
                entry.neighbour == kBoundaryStencil ? entry.boundaryValue : field[entry.neighbour];
            rhs += entry.weightedDelta * (neighbourValue - field[cell]);
        }
        gradients[cell] = gradientMatrixInverse_[cell].apply(rhs);
    }
    return gradients;
}

std::vector<Vector3> UnstructuredFvmBase::greenGaussGradients(const std::vector<double>& field) const {
    // grad(phi)_P = (1/V_P) * sum_f phi_f A_f_out. The face value is a plain
    // average for interior faces, the prescribed value where there is one,
    // and phi_P itself on a zero-gradient face (which is what zero gradient
    // means: the face carries the cell's own value).
    std::vector<Vector3> gradients(field.size());
    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        const auto& face = mesh_->face(f);
        if (mesh_->isBoundaryFace(f)) {
            double prescribed = 0.0;
            const bool hasValue = boundaryFaceValue_ && boundaryFaceValue_(f, prescribed);
            const double faceValue = hasValue ? prescribed : field[face.owner];
            gradients[face.owner] += face.areaVector * faceValue;
            continue;
        }
        const double faceValue = 0.5 * (field[face.owner] + field[face.neighbour]);
        gradients[face.owner] += face.areaVector * faceValue;
        gradients[face.neighbour] -= face.areaVector * faceValue;
    }
    for (std::size_t cell = 0; cell < gradients.size(); ++cell) {
        gradients[cell] = gradients[cell] / mesh_->cellVolume(cell);
    }
    return gradients;
}

std::vector<double> UnstructuredFvmBase::applyLaplacian(const std::vector<double>& x,
                                                         std::size_t pinnedCell) const {
    std::vector<double> result(x.size());
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = (cell == pinnedCell) ? x[cell] : laplacianDiagonal_[cell] * x[cell];
    }
    for (const InteriorFace& face : interiorFaces_) {
        if (face.owner != pinnedCell) {
            result[face.owner] -= face.coefficient * x[face.neighbour];
        }
        if (face.neighbour != pinnedCell) {
            result[face.neighbour] -= face.coefficient * x[face.owner];
        }
    }
    return result;
}

std::vector<double> UnstructuredFvmBase::nonOrthogonalCorrection(const std::vector<double>& field) const {
    // Evaluated from the previous iterate and moved to the right-hand side,
    // which is what keeps the matrix itself symmetric and positive definite.
    const std::vector<Vector3> gradients = computeCellGradients(field);

    std::vector<double> correction(field.size(), 0.0);
    for (const InteriorFace& face : interiorFaces_) {
        // **Corrected face gradient.** The plain average approximates the
        // gradient at the midpoint of the centroid segment, not at the face,
        // which is an O(h) error on a skewed mesh -- and since A_nonorth is
        // not perpendicular to d under the over-relaxed decomposition, that
        // error leaks straight into the flux. The fix is to keep the averaged
        // *tangential* part but replace the component along d by the compact
        // difference (phi_N - phi_P)/|d|, which is the one direction the two
        // cell values pin down accurately. Standard practice, and the
        // remaining first-order term identified when the least-squares
        // gradient alone left the order at ~1.
        const Vector3 averaged = averagedFaceGradient(face, gradients);
        const double compactNormalDerivative = (field[face.neighbour] - field[face.owner]) / face.distance;
        const Vector3 faceGradient =
            averaged + face.unitD * (compactNormalDerivative - averaged.dot(face.unitD));
        const double flux = faceGradient.dot(face.nonOrthogonalArea);
        correction[face.owner] += flux;
        correction[face.neighbour] -= flux;
    }
    for (const BoundaryFace& face : boundaryFaces_) {
        // Only the faces that entered the operator: a zero-gradient face
        // contributes no flux at all, orthogonal or not.
        if (!face.entersOperator) {
            continue;
        }
        double prescribed = 0.0;
        if (!boundaryFaceValue_ || !boundaryFaceValue_(face.meshFace, prescribed)) {
            continue;
        }
        const Vector3& cellGradient = gradients[face.cell];
        const double compactNormalDerivative = (prescribed - field[face.cell]) / face.distance;
        const Vector3 faceGradient =
            cellGradient + face.unitD * (compactNormalDerivative - cellGradient.dot(face.unitD));
        correction[face.cell] += faceGradient.dot(face.nonOrthogonalArea);
    }
    return correction;
}

std::size_t UnstructuredFvmBase::solveDeferredCorrection(std::vector<double>& x,
                                                          const std::vector<double>& baseRhs,
                                                          std::size_t pinnedCell,
                                                          std::size_t maxIterations, double tolerance,
                                                          std::size_t maxOuterSweeps,
                                                          double& lastOuterChange) const {
    const std::size_t n = x.size();
    std::size_t sweep = 0;
    for (; sweep < maxOuterSweeps; ++sweep) {
        // On the first sweep the correction comes from whatever `x` already
        // held -- zero for a fresh solve, which makes that sweep the
        // orthogonal-only solution, and the previous step's field for a
        // time-marching solver, which makes it consistent from the start.
        const std::vector<double> correction = nonOrthogonalCorrection(x);
        std::vector<double> rhs(n);
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] = baseRhs[i] + correction[i];
        }
        if (pinnedCell != kNoPinnedCell) {
            rhs[pinnedCell] = 0.0; // the pinned row is the identity; its value is fixed
        }

        const std::vector<double> previous = x;

        conjugateGradient(
            [this, pinnedCell](const std::vector<double>& v) { return applyLaplacian(v, pinnedCell); },
            rhs, x, maxIterations, tolerance);

        // Outer convergence: the correction has stopped moving the solution.
        double maxChange = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            maxChange = std::max(maxChange, std::fabs(x[i] - previous[i]));
        }
        lastOuterChange = maxChange;
        if (maxChange < tolerance) {
            break;
        }
    }
    return sweep;
}

double UnstructuredFvmBase::maxNonOrthogonality() const {
    double worst = 0.0;
    for (const InteriorFace& face : interiorFaces_) {
        worst = std::max(worst, face.nonOrthogonalArea.norm() / face.areaVector.norm());
    }
    return worst;
}

double UnstructuredFvmBase::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

} // namespace aether::solver
