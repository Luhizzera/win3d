#include "aether/solver/UnstructuredCavitySolver3D.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace aether::solver {

using core::Vector3;
using mesh::TetrahedralMesh;

UnstructuredCavitySolver3D::UnstructuredCavitySolver3D(
    const TetrahedralMesh& mesh, double viscosity,
    std::function<Vector3(const Vector3&)> wallVelocity,
    std::function<bool(const Vector3&)> isOutlet, double outletPressure)
    : UnstructuredFvmBase(mesh), viscosity_(viscosity), wallVelocity_(std::move(wallVelocity)),
      isOutlet_(std::move(isOutlet)), outletPressure_(outletPressure) {
    velocity_.assign(mesh.cellCount(), Vector3{});
    pressure_.assign(mesh.cellCount(), 0.0);

    // Only an outlet enters the Poisson operator. Pressure is prescribed
    // there, so unlike a wall the face is a real connection -- and it is what
    // makes the operator non-singular without pinning a reference cell. A
    // solid wall carries a zero-gradient pressure condition and therefore no
    // coefficient at all, exactly as the structured cavity solvers treat
    // their walls.
    buildFaceGeometry([this](std::size_t meshFace) {
        return isOutlet_ && isOutlet_(mesh_->face(meshFace).centroid);
    });
    buildBoundaryConditions();

    // Interior neighbours only. A solid wall carries a zero-gradient
    // condition for pressure, so including it as a stencil entry would
    // assert "the wall value equals the cell value", biasing the fit.
    buildGradientStencils();
    boundaryFlux_.assign(boundaryFaces_.size(), 0.0);
}

void UnstructuredCavitySolver3D::buildBoundaryConditions() {
    boundaryConditions_.assign(boundaryFaces_.size(), {});
    hasOutlet_ = false;
    maxWallSpeed_ = 0.0;
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        BoundaryCondition& condition = boundaryConditions_[i];
        condition.wallVelocity = wallVelocity_(boundaryFaces_[i].centroid);
        condition.isOutlet = isOutlet_ && isOutlet_(boundaryFaces_[i].centroid);
        maxWallSpeed_ = std::max(maxWallSpeed_, condition.wallVelocity.norm());
        hasOutlet_ = hasOutlet_ || condition.isOutlet;
    }
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
    const std::vector<Vector3> pressureGradients = leastSquaresGradients(pressure);
    std::vector<double> fluxes(interiorFaces_.size());
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const Vector3 interpolated = velocity[face.owner] * face.ownerWeight +
                                      velocity[face.neighbour] * (1.0 - face.ownerWeight);
        const Vector3 averagedGradient = averagedFaceGradient(face, pressureGradients);
        const double compactNormalDerivative =
            (pressure[face.neighbour] - pressure[face.owner]) / face.distance;
        const Vector3 correction =
            face.unitD * (averagedGradient.dot(face.unitD) - compactNormalDerivative);
        fluxes[i] = (interpolated + correction * dt).dot(face.areaVector);
    }
    return fluxes;
}

std::vector<double> UnstructuredCavitySolver3D::applyPoissonOperator(const std::vector<double>& x) const {
    // The same Laplacian Fase 2.2 validated -- minus its non-orthogonal
    // deferred correction, which this operator still lacks
    // (DIVIDA_TECNICA.md 1.2) -- with cell 0 pinned to remove the
    // pure-Neumann null space (every boundary here is a solid wall).
    // Cell 0 is pinned ONLY when every boundary is a wall. With an outlet
    // the prescribed pressure already fixes the level, and pinning as well
    // would over-constrain the system -- two incompatible references.
    const bool pinReference = !hasOutlet_;
    std::vector<double> result(x.size());
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = (pinReference && cell == 0) ? x[cell] : laplacianDiagonal_[cell] * x[cell];
    }
    for (const InteriorFace& face : interiorFaces_) {
        if (!pinReference || face.owner != 0) {
            result[face.owner] -= face.coefficient * x[face.neighbour];
        }
        if (!pinReference || face.neighbour != 0) {
            result[face.neighbour] -= face.coefficient * x[face.owner];
        }
    }
    return result;
}

std::vector<double> UnstructuredCavitySolver3D::applyHelmholtzOperator(const std::vector<double>& x,
                                                                        double dt) const {
    std::vector<double> result(x.size());
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = (mesh_->cellVolume(cell) / dt + viscosity_ * laplacianDiagonal_[cell]) * x[cell];
    }
    // A wall face stiffens its cell's equation even though the wall value
    // itself is known and belongs on the right-hand side.
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        if (boundaryConditions_[i].isOutlet) {
            continue; // zero-gradient: no viscous flux through an outlet
        }
        result[boundaryFaces_[i].cell] += viscosity_ * boundaryFaces_[i].coefficient * x[boundaryFaces_[i].cell];
    }
    for (const InteriorFace& face : interiorFaces_) {
        result[face.owner] -= viscosity_ * face.coefficient * x[face.neighbour];
        result[face.neighbour] -= viscosity_ * face.coefficient * x[face.owner];
    }
    return result;
}

std::vector<double> UnstructuredCavitySolver3D::solveHelmholtz(const std::vector<double>& rhs,
                                                                double dt) const {
    // Plain CG. This operator is far more diagonally dominant than the
    // pressure Poisson one -- the V/dt shift only adds to the diagonal -- so
    // it converges in very few iterations, which is why making diffusion
    // implicit costs much less than the step limit it removes.
    std::vector<double> x(rhs.size(), 0.0);
    conjugateGradient([this, dt](const std::vector<double>& v) { return applyHelmholtzOperator(v, dt); },
                      rhs, x, rhs.size(), 1e-12);
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
    // Boundary contributions: a prescribed-velocity face (wall or inlet)
    // injects its known flux, and an outlet both carries flux and pins the
    // pressure level through its Dirichlet coefficient.
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        const BoundaryFace& face = boundaryFaces_[i];
        const BoundaryCondition& condition = boundaryConditions_[i];
        const double flux = condition.isOutlet ? velocityStar[face.cell].dot(face.areaVector)
                                                : condition.wallVelocity.dot(face.areaVector);
        rhs[face.cell] -= flux / dt;
        if (condition.isOutlet) {
            rhs[face.cell] += face.coefficient * outletPressure_;
        }
    }
    if (!hasOutlet_) {
        rhs[0] = 0.0;
    }

    // Started from the previous step's pressure, which is why this converges
    // in a handful of iterations. With a pinned reference the operator's row
    // 0 is the identity and rhs[0] is zero, so pressure_[0] -- zero from
    // construction -- stays exactly zero and never enters any direction.
    conjugateGradient([this](const std::vector<double>& v) { return applyPoissonOperator(v); }, rhs,
                      pressure_, n, 1e-10);

    const std::vector<Vector3> pressureGradients = leastSquaresGradients(pressure_);
    for (std::size_t cell = 0; cell < n; ++cell) {
        velocity_[cell] = velocityStar[cell] - pressureGradients[cell] * dt;
    }

    // Record the boundary flux **the projection enforced**, using the same
    // compact face gradient the Poisson operator's outlet coefficient was
    // built from. Recomputing it later from the corrected cell velocity and
    // the least-squares gradient gives a different -- and systematically
    // wrong -- number, which is what produced the 13.2% imbalance.
    boundaryFlux_.assign(boundaryFaces_.size(), 0.0);
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        const BoundaryFace& face = boundaryFaces_[i];
        if (!boundaryConditions_[i].isOutlet) {
            continue;
        }
        boundaryFlux_[i] = velocityStar[face.cell].dot(face.areaVector) -
                            dt * face.coefficient * (outletPressure_ - pressure_[face.cell]);
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

    // Convective transport across the boundary. Without this, momentum enters
    // through an inlet and has no way to leave: it accumulates in the outlet
    // cells until the solution diverges. Found exactly that way -- the channel
    // case went to NaN while the closed cavity, which has no boundary flux at
    // all, was unaffected.
    for (std::size_t bi = 0; bi < boundaryFaces_.size(); ++bi) {
        const BoundaryFace& face = boundaryFaces_[bi];
        const BoundaryCondition& condition = boundaryConditions_[bi];
        const double massFlux =
            condition.isOutlet ? boundaryFlux_[bi] : condition.wallVelocity.dot(face.areaVector);
        // Upwind as everywhere else: outflow carries the cell value away,
        // inflow brings the prescribed boundary value in.
        const Vector3 upwind = massFlux >= 0.0 ? velocity_[face.cell] : condition.wallVelocity;
        flux[face.cell] -= upwind * massFlux;
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
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        const BoundaryCondition& condition = boundaryConditions_[i];
        if (condition.isOutlet) {
            continue;
        }
        const std::size_t cell = boundaryFaces_[i].cell;
        const double coefficient = viscosity_ * boundaryFaces_[i].coefficient;
        component[0][cell] += coefficient * condition.wallVelocity.x;
        component[1][cell] += coefficient * condition.wallVelocity.y;
        component[2][cell] += coefficient * condition.wallVelocity.z;
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

double UnstructuredCavitySolver3D::netBoundaryFlux() const {
    double net = 0.0;
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        // Inflow through a prescribed-velocity face counts too: that is how
        // an inlet enters the balance. Outlets use the flux the projection
        // actually enforced -- see boundaryFlux_.
        net += boundaryConditions_[i].isOutlet
                   ? boundaryFlux_[i]
                   : boundaryConditions_[i].wallVelocity.dot(boundaryFaces_[i].areaVector);
    }
    return net;
}

double UnstructuredCavitySolver3D::maxFaceDivergence() const {
    std::vector<double> divergence(velocity_.size(), 0.0);
    const std::vector<double> fluxes = faceMassFluxes(velocity_, pressure_, lastDt_);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        divergence[interiorFaces_[i].owner] += fluxes[i];
        divergence[interiorFaces_[i].neighbour] -= fluxes[i];
    }
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        divergence[boundaryFaces_[i].cell] +=
            boundaryConditions_[i].isOutlet
                ? boundaryFlux_[i]
                : boundaryConditions_[i].wallVelocity.dot(boundaryFaces_[i].areaVector);
    }
    double worst = 0.0;
    for (std::size_t cell = 0; cell < divergence.size(); ++cell) {
        worst = std::max(worst, std::fabs(divergence[cell]) / mesh_->cellVolume(cell));
    }
    return worst;
}

} // namespace aether::solver
