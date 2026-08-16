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
    // **Convection only.** With the viscous term solved implicitly there is
    // no diffusive stability bound at all, which is what that change was
    // for: the explicit bound scaled with the square of the smallest cell,
    // and a Delaunay tetrahedralization always produces slivers, so the
    // mesh's worst cell -- not the physics -- was setting the step.
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
    // Safety factor, deliberately: ROADMAP Fase 1 found the structured cavity
    // running at CFL exactly 1.0000 with no margin at all, which is why a
    // small perturbation tipped it into NaN.
    return 0.4 * smallestScale / std::max(maxWallSpeed_, 1e-12);
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

std::vector<double> UnstructuredCavitySolver3D::applyHelmholtzOperator(const std::vector<double>& x,
                                                                        double dt) const {
    std::vector<double> result(x.size());
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = (mesh_->cellVolume(cell) / dt + viscosity_ * poissonDiagonal_[cell]) * x[cell];
    }
    // A wall face stiffens its cell's equation even though the wall value
    // itself is known and belongs on the right-hand side.
    for (const BoundaryFace& face : boundaryFaces_) {
        result[face.cell] += viscosity_ * face.laplacianCoefficient * x[face.cell];
    }
    for (const InteriorFace& face : interiorFaces_) {
        result[face.owner] -= viscosity_ * face.laplacianCoefficient * x[face.neighbour];
        result[face.neighbour] -= viscosity_ * face.laplacianCoefficient * x[face.owner];
    }
    return result;
}

std::vector<double> UnstructuredCavitySolver3D::solveHelmholtz(const std::vector<double>& rhs,
                                                                double dt) const {
    // Plain CG. This operator is far more diagonally dominant than the
    // pressure Poisson one -- the V/dt shift only adds to the diagonal -- so
    // it converges in very few iterations, which is why making diffusion
    // implicit costs much less than the step limit it removes.
    const std::size_t n = rhs.size();
    std::vector<double> x(n, 0.0);
    std::vector<double> residual = rhs;
    std::vector<double> direction = residual;
    double residualDotResidual = 0.0;
    for (double r : residual) {
        residualDotResidual += r * r;
    }
    for (std::size_t iteration = 0; iteration < n; ++iteration) {
        if (std::sqrt(residualDotResidual) < 1e-12) {
            break;
        }
        const std::vector<double> applied = applyHelmholtzOperator(direction, dt);
        double directionDotApplied = 0.0;
        for (std::size_t i2 = 0; i2 < n; ++i2) {
            directionDotApplied += direction[i2] * applied[i2];
        }
        if (directionDotApplied == 0.0) {
            break;
        }
        const double alpha = residualDotResidual / directionDotApplied;
        double newResidualDotResidual = 0.0;
        for (std::size_t i2 = 0; i2 < n; ++i2) {
            x[i2] += alpha * direction[i2];
            residual[i2] -= alpha * applied[i2];
            newResidualDotResidual += residual[i2] * residual[i2];
        }
        const double beta = newResidualDotResidual / residualDotResidual;
        for (std::size_t i2 = 0; i2 < n; ++i2) {
            direction[i2] = residual[i2] + beta * direction[i2];
        }
        residualDotResidual = newResidualDotResidual;
    }
    return x;
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

    // --- Viscous diffusion, **implicit**. Per component,
    //   (V/dt + nu L) u* = (V/dt) (u^n + dt * convection / V) + nu * wall terms
    // so no diffusive stability bound applies at all. The three components
    // are independent scalar solves: the viscous operator does not couple
    // them.
    std::vector<std::vector<double>> component(3, std::vector<double>(n));
    for (std::size_t cell = 0; cell < n; ++cell) {
        const double volumeOverDt = mesh_->cellVolume(cell) / dt;
        const Vector3 explicitPart = velocity_[cell] + flux[cell] * (dt / mesh_->cellVolume(cell));
        component[0][cell] = volumeOverDt * explicitPart.x;
        component[1][cell] = volumeOverDt * explicitPart.y;
        component[2][cell] = volumeOverDt * explicitPart.z;
    }
    // The wall velocity is known, so its viscous flux is a source rather than
    // part of the operator. This is where the lid drives the flow.
    for (const BoundaryFace& face : boundaryFaces_) {
        const double coefficient = viscosity_ * face.laplacianCoefficient;
        component[0][face.cell] += coefficient * face.wallVelocity.x;
        component[1][face.cell] += coefficient * face.wallVelocity.y;
        component[2][face.cell] += coefficient * face.wallVelocity.z;
    }

    std::vector<Vector3> velocityStar(n);
    for (int c = 0; c < 3; ++c) {
        const std::vector<double> solved = solveHelmholtz(component[c], dt);
        for (std::size_t cell = 0; cell < n; ++cell) {
            if (c == 0) {
                velocityStar[cell].x = solved[cell];
            } else if (c == 1) {
                velocityStar[cell].y = solved[cell];
            } else {
                velocityStar[cell].z = solved[cell];
            }
        }
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
