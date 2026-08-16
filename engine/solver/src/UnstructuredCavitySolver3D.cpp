#include "aether/solver/UnstructuredCavitySolver3D.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace aether::solver {

using core::Vector3;
using mesh::TetrahedralMesh;

UnstructuredCavitySolver3D::UnstructuredCavitySolver3D(
    const TetrahedralMesh& mesh, double viscosity,
    std::function<Vector3(const Vector3&)> wallVelocity)
    : mesh_(&mesh), viscosity_(viscosity), wallVelocity_(std::move(wallVelocity)) {
    velocity_.assign(mesh.cellCount(), Vector3{});
    pressure_.assign(mesh.cellCount(), 0.0);
    buildFaces();
    buildGradientStencils();
}

void UnstructuredCavitySolver3D::buildFaces() {
    poissonDiagonal_.assign(mesh_->cellCount(), 0.0);

    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        const auto& face = mesh_->face(f);

        if (mesh_->isBoundaryFace(f)) {
            const Vector3 d = face.centroid - mesh_->cellCentroid(face.owner);
            const double areaDotD = face.areaVector.dot(d);
            if (areaDotD <= 0.0) {
                continue;
            }
            BoundaryFace boundary;
            boundary.cell = face.owner;
            boundary.areaVector = face.areaVector;
            boundary.centroid = face.centroid;
            boundary.wallVelocity = wallVelocity_(face.centroid);
            boundary.laplacianCoefficient = face.areaVector.normSquared() / areaDotD;
            boundary.distance = d.norm();
            maxWallSpeed_ = std::max(maxWallSpeed_, boundary.wallVelocity.norm());
            boundaryFaces_.push_back(boundary);
            // No pressure coefficient: the pressure boundary condition is
            // zero-gradient at a solid wall, so a boundary face contributes
            // nothing to the Poisson operator -- exactly as the structured
            // cavity solvers treat their walls.
            continue;
        }

        const Vector3 d = mesh_->cellCentroid(face.neighbour) - mesh_->cellCentroid(face.owner);
        const double areaDotD = face.areaVector.dot(d);
        if (areaDotD <= 0.0) {
            continue;
        }

        InteriorFace interior;
        interior.owner = face.owner;
        interior.neighbour = face.neighbour;
        interior.laplacianCoefficient = face.areaVector.normSquared() / areaDotD;
        interior.areaVector = face.areaVector;
        interior.nonOrthogonalArea = face.areaVector - d * interior.laplacianCoefficient;
        interior.distance = d.norm();
        interior.unitD = d / interior.distance;

        const Vector3 unitNormal = face.areaVector / face.areaVector.norm();
        const double ownerDistance =
            std::fabs((face.centroid - mesh_->cellCentroid(face.owner)).dot(unitNormal));
        const double neighbourDistance =
            std::fabs((mesh_->cellCentroid(face.neighbour) - face.centroid).dot(unitNormal));
        const double sum = ownerDistance + neighbourDistance;
        interior.ownerWeight = sum > 0.0 ? neighbourDistance / sum : 0.5;

        poissonDiagonal_[face.owner] += interior.laplacianCoefficient;
        poissonDiagonal_[face.neighbour] += interior.laplacianCoefficient;
        interiorFaces_.push_back(interior);
    }
}

void UnstructuredCavitySolver3D::buildGradientStencils() {
    const std::size_t n = mesh_->cellCount();
    gradientStencil_.assign(n, {});
    gradientMatrixInverse_.assign(n, {});

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

    // Interior neighbours only. A solid wall carries a zero-gradient
    // condition for pressure, so including it as a stencil entry would
    // assert "the wall value equals the cell value", biasing the fit.
    for (const InteriorFace& face : interiorFaces_) {
        const Vector3 d = face.unitD * face.distance;
        const double w = 1.0 / (face.distance * face.distance);
        gradientStencil_[face.owner].push_back({face.neighbour, d * w});
        gradientStencil_[face.neighbour].push_back({face.owner, -d * w});
        accumulate(matrices[face.owner], d, w);
        accumulate(matrices[face.neighbour], d, w);
    }

    for (std::size_t cell = 0; cell < n; ++cell) {
        const Sym& m = matrices[cell];
        const double c11 = m.yy * m.zz - m.yz * m.yz;
        const double c12 = m.xz * m.yz - m.xy * m.zz;
        const double c13 = m.xy * m.yz - m.yy * m.xz;
        const double det = m.xx * c11 + m.xy * c12 + m.xz * c13;
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

std::vector<Vector3> UnstructuredCavitySolver3D::scalarGradients(const std::vector<double>& field) const {
    std::vector<Vector3> gradients(field.size());
    for (std::size_t cell = 0; cell < field.size(); ++cell) {
        if (!gradientMatrixInverse_[cell].valid) {
            continue; // leaves a zero gradient, the safe choice for a degenerate stencil
        }
        Vector3 rhs;
        for (const GradientStencilEntry& entry : gradientStencil_[cell]) {
            rhs += entry.weightedDelta * (field[entry.neighbour] - field[cell]);
        }
        gradients[cell] = gradientMatrixInverse_[cell].apply(rhs);
    }
    return gradients;
}

double UnstructuredCavitySolver3D::stableTimeStep() const {
    // **The limit is taken from the operator itself, not from a length-scale
    // proxy.** The explicit viscous update is
    //     u_P += dt * nu * sum_f a_f (u_N - u_P) / V_P
    // whose diagonal coefficient is nu * sum_f a_f / V_P, so stability needs
    // dt below V_P / (nu * sum_f a_f) for every cell. That is exactly
    // computable here because the a_f are already assembled.
    //
    // A first version instead used (V_P / total face area)^2 / (6 nu) as a
    // stand-in for the cell's size. It was badly pessimistic on the slivers
    // every Delaunay tetrahedralization produces -- volume shrinks toward
    // zero while face area does not, so the proxy collapses much faster than
    // the true limit does. Measured directly: that made a single cavity test
    // take 6m57s. The proxy was never the criterion, only a guess at it.
    double diffusiveLimit = 1e300;
    if (viscosity_ > 0.0) {
        std::vector<double> coefficientSum(mesh_->cellCount(), 0.0);
        for (const InteriorFace& face : interiorFaces_) {
            coefficientSum[face.owner] += face.laplacianCoefficient;
            coefficientSum[face.neighbour] += face.laplacianCoefficient;
        }
        for (const BoundaryFace& face : boundaryFaces_) {
            coefficientSum[face.cell] += face.laplacianCoefficient;
        }
        for (std::size_t cell = 0; cell < mesh_->cellCount(); ++cell) {
            if (coefficientSum[cell] > 0.0) {
                diffusiveLimit =
                    std::min(diffusiveLimit, mesh_->cellVolume(cell) / (viscosity_ * coefficientSum[cell]));
            }
        }
    }

    // Convection still needs a length scale, and here volume / face area is
    // the right one: it is the distance the flux actually has to cross.
    double smallestScale = 1e300;
    for (std::size_t cell = 0; cell < mesh_->cellCount(); ++cell) {
        double areaSum = 0.0;
        for (std::size_t faceIndex : mesh_->cellFaces(cell)) {
            areaSum += mesh_->face(faceIndex).areaVector.norm();
        }
        if (areaSum > 0.0) {
            smallestScale = std::min(smallestScale, mesh_->cellVolume(cell) / areaSum);
        }
    }
    const double convectiveLimit = smallestScale / std::max(maxWallSpeed_, 1e-12);

    // Safety factor, deliberately: ROADMAP Fase 1 found the structured
    // cavity running at CFL exactly 1.0000 with no margin at all, which is
    // why a small perturbation tipped it into NaN.
    //
    // **This factor is also the knob on a measured tradeoff.** Switching from
    // the length-scale proxy to the criterion above cut the cavity test from
    // 58.8s to 21.8s with the flow topology unchanged (u mean +0.0439/-0.0099
    // before, +0.0457/-0.0098 after), but face divergence rose from 3.3e-04
    // to 6.4e-02: a larger step leaves more for each projection to remove.
    // Both are honest consequences of a bigger dt, not a defect introduced by
    // the criterion. Halve this factor if mass conservation matters more than
    // wall-clock for a given run -- it stays far cheaper than the proxy was.
    return 0.4 * std::min(diffusiveLimit, convectiveLimit);
}

std::vector<double> UnstructuredCavitySolver3D::faceMassFluxes(const std::vector<Vector3>& velocity,
                                                                 const std::vector<double>& pressure,
                                                                 double dt) const {
    // Rhie-Chow style: interpolate the velocity to the face, then replace the
    // pressure-gradient part of that interpolation with the compact face
    // difference. Without this the face flux is blind to a cell-to-cell
    // pressure oscillation, which is the collocated-grid failure mode this
    // arrangement is otherwise exposed to.
    const std::vector<Vector3> pressureGradients = scalarGradients(pressure);
    std::vector<double> fluxes(interiorFaces_.size());
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const Vector3 interpolated = velocity[face.owner] * face.ownerWeight +
                                      velocity[face.neighbour] * (1.0 - face.ownerWeight);
        const Vector3 averagedGradient = pressureGradients[face.owner] * face.ownerWeight +
                                          pressureGradients[face.neighbour] * (1.0 - face.ownerWeight);
        const double compactNormalDerivative =
            (pressure[face.neighbour] - pressure[face.owner]) / face.distance;
        const Vector3 correction =
            face.unitD * (averagedGradient.dot(face.unitD) - compactNormalDerivative);
        fluxes[i] = (interpolated + correction * dt).dot(face.areaVector);
    }
    return fluxes;
}

std::vector<double> UnstructuredCavitySolver3D::applyPoissonOperator(const std::vector<double>& x) const {
    // Same operator Fase 2.2 validated, with cell 0 pinned to remove the
    // pure-Neumann null space (every boundary here is a solid wall).
    std::vector<double> result(x.size());
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = cell == 0 ? x[cell] : poissonDiagonal_[cell] * x[cell];
    }
    for (const InteriorFace& face : interiorFaces_) {
        if (face.owner != 0) {
            result[face.owner] -= face.laplacianCoefficient * x[face.neighbour];
        }
        if (face.neighbour != 0) {
            result[face.neighbour] -= face.laplacianCoefficient * x[face.owner];
        }
    }
    return result;
}

void UnstructuredCavitySolver3D::projectToDivergenceFree(std::vector<Vector3>& velocityStar, double dt) {
    const std::size_t n = velocityStar.size();

    // Right-hand side: the face-flux divergence of the predicted velocity.
    // Plain interpolation here (dt = 0 in the Rhie-Chow term) because the
    // predictor carries no pressure correction to be consistent with.
    const std::vector<double> starFluxes = faceMassFluxes(velocityStar, pressure_, 0.0);
    std::vector<double> rhs(n, 0.0);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        rhs[interiorFaces_[i].owner] -= starFluxes[i] / dt;
        rhs[interiorFaces_[i].neighbour] += starFluxes[i] / dt;
    }
    // Solid walls: no flow penetrates, so a boundary face contributes no
    // mass flux at all.
    rhs[0] = 0.0;

    std::vector<double> residual(n);
    {
        const std::vector<double> applied = applyPoissonOperator(pressure_);
        for (std::size_t i = 0; i < n; ++i) {
            residual[i] = i == 0 ? 0.0 : rhs[i] - applied[i];
        }
    }
    std::vector<double> direction = residual;
    double residualDotResidual = 0.0;
    for (double r : residual) {
        residualDotResidual += r * r;
    }

    for (std::size_t iteration = 0; iteration < n; ++iteration) {
        if (std::sqrt(residualDotResidual) < 1e-10) {
            break;
        }
        const std::vector<double> applied = applyPoissonOperator(direction);
        double directionDotApplied = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            directionDotApplied += direction[i] * applied[i];
        }
        if (directionDotApplied == 0.0) {
            break;
        }
        const double alpha = residualDotResidual / directionDotApplied;
        double newResidualDotResidual = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (i != 0) {
                pressure_[i] += alpha * direction[i];
                residual[i] -= alpha * applied[i];
            }
            newResidualDotResidual += residual[i] * residual[i];
        }
        const double beta = newResidualDotResidual / residualDotResidual;
        for (std::size_t i = 0; i < n; ++i) {
            direction[i] = i == 0 ? 0.0 : residual[i] + beta * direction[i];
        }
        residualDotResidual = newResidualDotResidual;
    }

    const std::vector<Vector3> pressureGradients = scalarGradients(pressure_);
    for (std::size_t cell = 0; cell < n; ++cell) {
        velocity_[cell] = velocityStar[cell] - pressureGradients[cell] * dt;
    }
}

void UnstructuredCavitySolver3D::step(double dt) {
    const std::size_t n = velocity_.size();
    std::vector<Vector3> flux(n, Vector3{});

    // --- Convection, upwinded on the face mass flux (see the class comment
    // for why not central differencing), in conservative form so that what
    // leaves one cell enters its neighbour exactly.
    const std::vector<double> massFluxes = faceMassFluxes(velocity_, pressure_, lastDt_);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const double massFlux = massFluxes[i];
        const Vector3 upwind = massFlux >= 0.0 ? velocity_[face.owner] : velocity_[face.neighbour];
        flux[face.owner] -= upwind * massFlux;
        flux[face.neighbour] += upwind * massFlux;
    }

    // --- Viscous diffusion, orthogonal part plus the non-orthogonal
    // correction evaluated from the current field (same decomposition as
    // UnstructuredDiffusionSolver).
    for (const InteriorFace& face : interiorFaces_) {
        const Vector3 difference = velocity_[face.neighbour] - velocity_[face.owner];
        const Vector3 viscous = difference * (viscosity_ * face.laplacianCoefficient);
        flux[face.owner] += viscous;
        flux[face.neighbour] -= viscous;
    }
    // Wall viscous flux: the prescribed wall velocity acts as the "neighbour"
    // value at the face centroid. This is what drives the whole flow -- the
    // lid's tangential motion enters the momentum equation here and nowhere
    // else.
    for (const BoundaryFace& face : boundaryFaces_) {
        const Vector3 difference = face.wallVelocity - velocity_[face.cell];
        flux[face.cell] += difference * (viscosity_ * face.laplacianCoefficient);
    }

    std::vector<Vector3> velocityStar(n);
    for (std::size_t cell = 0; cell < n; ++cell) {
        velocityStar[cell] = velocity_[cell] + flux[cell] * (dt / mesh_->cellVolume(cell));
    }

    projectToDivergenceFree(velocityStar, dt);
    lastDt_ = dt;
    time_ += dt;
}

double UnstructuredCavitySolver3D::maxFaceDivergence() const {
    std::vector<double> divergence(velocity_.size(), 0.0);
    const std::vector<double> fluxes = faceMassFluxes(velocity_, pressure_, lastDt_);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        divergence[interiorFaces_[i].owner] += fluxes[i];
        divergence[interiorFaces_[i].neighbour] -= fluxes[i];
    }
    double worst = 0.0;
    for (std::size_t cell = 0; cell < divergence.size(); ++cell) {
        worst = std::max(worst, std::fabs(divergence[cell]) / mesh_->cellVolume(cell));
    }
    return worst;
}

} // namespace aether::solver
