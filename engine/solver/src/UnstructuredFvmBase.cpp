#include "aether/solver/UnstructuredFvmBase.hpp"

#include "aether/solver/ConvectionLimiter.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

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

    refreshBoundaryValues();
}

void UnstructuredFvmBase::refreshBoundaryValues() {
    boundaryHasValue_.assign(mesh_->faceCount(), 0);
    boundaryValueCache_.assign(mesh_->faceCount(), 0.0);
    if (!boundaryFaceValue_) {
        return;
    }
    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        if (!mesh_->isBoundaryFace(f)) {
            continue;
        }
        double value = 0.0;
        if (boundaryFaceValue_(f, value)) {
            boundaryHasValue_[f] = 1;
            boundaryValueCache_[f] = value;
        }
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
            if (!useBoundaryValues || boundaryHasValue_[f] == 0) {
                continue; // zero-gradient face: carries no gradient information
            }
            const double prescribed = boundaryValueCache_[f];
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
        gradientMatrixInverse_[cell] = invertSymmetric3x3(m.xx, m.xy, m.xz, m.yy, m.yz, m.zz);
        if (!gradientMatrixInverse_[cell].valid) {
            ++deficientStencilCount_;
        }
    }
}

UnstructuredFvmBase::SymmetricInverse UnstructuredFvmBase::invertSymmetric3x3(double xx, double xy,
                                                                              double xz, double yy,
                                                                              double yz, double zz) {
    // Cofactors of a symmetric 3x3.
    const double c11 = yy * zz - yz * yz;
    const double c12 = xz * yz - xy * zz;
    const double c13 = xy * yz - yy * xz;
    const double det = xx * c11 + xy * c12 + xz * c13;

    SymmetricInverse inverse;
    const double scale = xx + yy + zz;
    if (scale <= 0.0 || std::fabs(det) < 1e-12 * scale * scale * scale) {
        inverse.valid = false;
        return inverse;
    }

    const double inverseDet = 1.0 / det;
    inverse.xx = c11 * inverseDet;
    inverse.xy = c12 * inverseDet;
    inverse.xz = c13 * inverseDet;
    inverse.yy = (xx * zz - xz * xz) * inverseDet;
    inverse.yz = (xy * xz - xx * yz) * inverseDet;
    inverse.zz = (xx * yy - xy * xy) * inverseDet;
    inverse.valid = true;
    return inverse;
}

std::vector<Vector3> UnstructuredFvmBase::computeCellGradients(const std::vector<double>& field) const {
    // The Green-Gauss fallback is only built when this mesh actually has a
    // cell that needs it -- it costs two full passes over the faces.
    std::vector<Vector3> fallback;
    std::vector<double> steepestSlope;
    if (deficientStencilCount_ > 0) {
        fallback = greenGaussGradients(field);
        // **The bound the fallback is held to, and why it has to exist.**
        // Green-Gauss on a cell with one or two interior neighbours is not a
        // gradient estimate so much as a guess: its magnitude scales with
        // |A|/V and nothing constrains it. That was measured the hard way --
        // replacing the old silent zero with Green-Gauss (DIVIDA_TECNICA.md
        // 1.3) is what pushed the spectral radius of a single step from
        // 0.9776 to 2.0611 on a distorted mesh, which is the whole of item
        // 4.3. The instability lived on exactly these cells: 88.8% of the
        // unstable mode's energy sat on fallback cells, against 3-7% on
        // meshes that run.
        //
        // The bound is the steepest slope the cell can actually see: no
        // neighbour difference divided by its distance. A gradient larger
        // than every difference it was built from is not extrapolation, it is
        // invention. Where the mesh is fine the bound never binds and nothing
        // changes; where it binds, the estimate is pulled back to the largest
        // rate the data supports, which keeps it a bounded operator.
        steepestSlope.assign(field.size(), 0.0);
        for (const InteriorFace& face : interiorFaces_) {
            const double slope = std::fabs(field[face.neighbour] - field[face.owner]) / face.distance;
            steepestSlope[face.owner] = std::max(steepestSlope[face.owner], slope);
            steepestSlope[face.neighbour] = std::max(steepestSlope[face.neighbour], slope);
        }
        for (const BoundaryFace& face : boundaryFaces_) {
            if (boundaryHasValue_[face.meshFace] == 0) {
                continue; // a zero-gradient face supports no slope at all
            }
            const double slope =
                std::fabs(boundaryValueCache_[face.meshFace] - field[face.cell]) / face.distance;
            steepestSlope[face.cell] = std::max(steepestSlope[face.cell], slope);
        }
    }

    std::vector<Vector3> gradients(field.size());
    for (std::size_t cell = 0; cell < field.size(); ++cell) {
        if (!gradientMatrixInverse_[cell].valid) {
            const Vector3& candidate = fallback[cell];
            const double magnitude = candidate.norm();
            gradients[cell] = magnitude > steepestSlope[cell] && magnitude > 0.0
                                  ? candidate * (steepestSlope[cell] / magnitude)
                                  : candidate;
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
            const double faceValue =
                boundaryHasValue_[f] != 0 ? boundaryValueCache_[f] : field[face.owner];
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

double UnstructuredFvmBase::faceValue(const InteriorFace& face, double massFlux,
                                       const std::vector<double>& field,
                                       const std::vector<Vector3>& gradients,
                                       ConvectionScheme scheme) const {
    const bool ownerIsUpwind = massFlux >= 0.0;
    const std::size_t upwind = ownerIsUpwind ? face.owner : face.neighbour;
    if (scheme == ConvectionScheme::FirstOrderUpwind) {
        return field[upwind];
    }

    const std::size_t downwind = ownerIsUpwind ? face.neighbour : face.owner;
    const double difference = field[downwind] - field[upwind];

    // Distance-weighted interpolation, the second-order end of the blend.
    const double central =
        field[face.owner] * face.ownerWeight + field[face.neighbour] * (1.0 - face.ownerWeight);

    if (faceDifferenceIsNegligible(difference, field[upwind], field[downwind])) {
        return field[upwind];
    }

    // d from the upwind cell to the downwind one, so the ratio is measured
    // along the direction the flux actually travels.
    const Vector3 upwindToDownwind = ownerIsUpwind ? face.unitD * face.distance
                                                    : face.unitD * (-face.distance);
    const double ratio = 2.0 * gradients[upwind].dot(upwindToDownwind) / difference - 1.0;
    return field[upwind] + vanLeerLimiter(ratio) * (central - field[upwind]);
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
        if (boundaryHasValue_[face.meshFace] == 0) {
            continue;
        }
        const double prescribed = boundaryValueCache_[face.meshFace];
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
                                                          double& lastOuterChange,
                                                          double& relaxation) const {
    const std::size_t n = x.size();
    // **The deferred correction is a fixed-point iteration, and it is not
    // always a contraction.** On a mesh with non-orthogonality around 15 it
    // diverges outright: the per-sweep change reaches 7e+07 and the returned
    // field 1e+08, which this function used to hand back as a solution
    // (DIVIDA_TECNICA.md 3.3). That is the same failure the Navier-Stokes
    // step already refuses -- a wrong answer that looks like a result -- and
    // it is worse here, because a caller reading a convergence order off a
    // diverging solve gets a number that means nothing at all.
    //
    // Detected by comparing against the best sweep so far rather than the
    // first: the correction legitimately makes things worse for a sweep or
    // two before settling, so growth from the *start* is not evidence, while
    // growth of three orders of magnitude away from the best point ever
    // reached is not something a converging iteration does.
    double bestChange = std::numeric_limits<double>::infinity();
    // **Under-relaxation, applied only when the iteration misbehaves.**
    // Deferred correction is a fixed point whose contraction factor is
    // ρ(A⁻¹C); on a bad enough mesh that exceeds 1 and the iteration diverges
    // -- measured at three separate places in this project, and no static mesh
    // metric predicted any of them (non-orthogonality, the non-orthogonal to
    // orthogonal area ratio, and the deficient-stencil count were all tried
    // and all failed to separate a mesh that runs from one that does not).
    //
    // So the damping is adaptive rather than mesh-derived: applying α to the
    // correction scales the contraction factor by α, so halving until the
    // iteration behaves finds a working α without needing to predict it. It
    // engages only after the change has grown ten-fold past its best value,
    // which a converging iteration never does -- the correction legitimately
    // makes things worse for a sweep or two, and damping that would be
    // slowing down a solve that was fine.
    constexpr double kMinimumRelaxation = 1.0 / 64.0;
    std::size_t sweep = 0;
    for (; sweep < maxOuterSweeps; ++sweep) {
        // On the first sweep the correction comes from whatever `x` already
        // held -- zero for a fresh solve, which makes that sweep the
        // orthogonal-only solution, and the previous step's field for a
        // time-marching solver, which makes it consistent from the start.
        const std::vector<double> correction = nonOrthogonalCorrection(x);
        std::vector<double> rhs(n);
        for (std::size_t i = 0; i < n; ++i) {
            rhs[i] = baseRhs[i] + relaxation * correction[i];
        }
        if (pinnedCell != kNoPinnedCell) {
            rhs[pinnedCell] = 0.0; // the pinned row is the identity; its value is fixed
        }

        const std::vector<double> previous = x;

        conjugateGradient(
            [this, pinnedCell](const std::vector<double>& v) { return applyLaplacian(v, pinnedCell); },
            rhs, x, maxIterations, tolerance);

        // Outer convergence: the correction has stopped moving the solution.
        // **Relative, not absolute, and that distinction was measured.**
        // `tolerance` used to be compared against the raw change, so a
        // pressure field of order 100 needed to settle to 1e-10 in absolute
        // terms -- a criterion that essentially never fired, leaving every
        // solve to run its full sweep cap (DIVIDA_TECNICA.md 5.4). Scaling by
        // the field's own magnitude makes the same number mean the same thing
        // whatever the units are, which is what a caller passing 1e-10
        // expects it to mean.
        //
        // The absolute form is kept as a floor for the case the relative one
        // cannot express: a field that is identically zero has no scale, and
        // a rest state is a legitimate solution.
        double maxChange = 0.0;
        double scale = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            maxChange = std::max(maxChange, std::fabs(x[i] - previous[i]));
            scale = std::max(scale, std::fabs(x[i]));
        }
        lastOuterChange = maxChange;
        if (maxChange < tolerance || maxChange < tolerance * scale) {
            break;
        }
        bestChange = std::min(bestChange, maxChange);
        if (!std::isfinite(maxChange) || maxChange > 10.0 * bestChange) {
            if (relaxation > kMinimumRelaxation) {
                // Diverging at this α: damp and let it try again. `bestChange`
                // is reset because the iteration being measured is now a
                // different one, and judging the new α against the old one's
                // best would trip the guard immediately.
                relaxation *= 0.5;
                bestChange = std::numeric_limits<double>::infinity();
                continue;
            }
            throw std::runtime_error(
                "UnstructuredFvmBase: the non-orthogonal deferred correction is diverging even "
                "under-relaxed to 1/64 -- the per-sweep change keeps growing past ten times its "
                "best value at every damping level tried. This is a property of the mesh, not of "
                "the sweep cap: raising the cap makes it worse. See DIVIDA_TECNICA.md 3.3.");
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
