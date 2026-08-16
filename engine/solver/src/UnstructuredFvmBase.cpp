#include "aether/solver/UnstructuredFvmBase.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

using core::Vector3;

void UnstructuredFvmBase::buildFaceGeometry(
    const std::function<bool(std::size_t)>& boundaryEntersOperator) {
    interiorFaces_.clear();
    boundaryFaces_.clear();
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
            boundaryFaces_.push_back(boundary);
            if (boundaryEntersOperator && boundaryEntersOperator(f)) {
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

        laplacianDiagonal_[face.owner] += coefficient;
        laplacianDiagonal_[face.neighbour] += coefficient;
        interiorFaces_.push_back(interior);
    }
}

void UnstructuredFvmBase::buildGradientStencils(
    const std::function<bool(std::size_t, double&)>& prescribedBoundaryValue) {
    const std::size_t n = mesh_->cellCount();
    gradientStencil_.assign(n, {});
    gradientMatrixInverse_.assign(n, {});
    hasDeficientStencil_ = false;

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
            if (!prescribedBoundaryValue || !prescribedBoundaryValue(f, prescribed)) {
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
            hasDeficientStencil_ = true;
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

std::vector<Vector3> UnstructuredFvmBase::leastSquaresGradients(
    const std::vector<double>& field, const std::vector<Vector3>* fallback) const {
    std::vector<Vector3> gradients(field.size());
    for (std::size_t cell = 0; cell < field.size(); ++cell) {
        if (!gradientMatrixInverse_[cell].valid) {
            if (fallback != nullptr) {
                gradients[cell] = (*fallback)[cell];
            }
            continue; // otherwise left at zero -- see the header on why that is a choice
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
