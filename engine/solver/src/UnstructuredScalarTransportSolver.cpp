#include "aether/solver/UnstructuredScalarTransportSolver.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

using core::Vector3;
using mesh::TetrahedralMesh;

UnstructuredScalarTransportSolver::UnstructuredScalarTransportSolver(
    const TetrahedralMesh& mesh, double diffusivity,
    std::function<Vector3(const Vector3&)> velocity, ConvectionScheme scheme)
    : UnstructuredFvmBase(mesh), diffusivity_(diffusivity), velocity_(std::move(velocity)),
      scheme_(scheme) {
    phi_.assign(mesh.cellCount(), 0.0);
    boundaryIsDirichlet_.assign(mesh.faceCount(), false);
    boundaryValue_.assign(mesh.faceCount(), 0.0);
    setBoundaryFaceValue([this](std::size_t meshFace, double& value) {
        if (!boundaryIsDirichlet_[meshFace]) {
            return false;
        }
        value = boundaryValue_[meshFace];
        return true;
    });
    buildFaceGeometry([this](std::size_t meshFace) { return boundaryIsDirichlet_[meshFace]; });
    buildGradientStencils(/*useBoundaryValues=*/true);

    // The velocity is prescribed, so its face fluxes are geometry: built once.
    interiorMassFlux_.assign(interiorFaces_.size(), 0.0);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        // Sampled at the point where the centroid segment crosses the face
        // plane, weighted the same way every other face quantity here is,
        // rather than at the face centroid: the two differ on a skewed mesh
        // and this is the point the face's flux is attributed to.
        const Vector3 point = mesh_->cellCentroid(face.owner) +
                              face.unitD * (face.distance * (1.0 - face.ownerWeight));
        interiorMassFlux_[i] = velocity_(point).dot(face.areaVector);
    }
    boundaryMassFlux_.assign(boundaryFaces_.size(), 0.0);
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        boundaryMassFlux_[i] = velocity_(boundaryFaces_[i].centroid).dot(boundaryFaces_[i].areaVector);
    }
}

void UnstructuredScalarTransportSolver::setDirichletBoundary(
    const std::function<bool(const Vector3&)>& selector,
    const std::function<double(const Vector3&)>& value) {
    for (std::size_t f = 0; f < mesh_->faceCount(); ++f) {
        if (!mesh_->isBoundaryFace(f)) {
            continue;
        }
        const Vector3& centroid = mesh_->face(f).centroid;
        if (selector(centroid)) {
            boundaryIsDirichlet_[f] = true;
            boundaryValue_[f] = value(centroid);
        }
    }
    buildFaceGeometry([this](std::size_t meshFace) { return boundaryIsDirichlet_[meshFace]; });
    buildGradientStencils(/*useBoundaryValues=*/true);
}

void UnstructuredScalarTransportSolver::setSourceTerm(
    const std::function<double(const Vector3&)>& source) {
    integratedSource_.assign(mesh_->cellCount(), 0.0);
    for (std::size_t cell = 0; cell < mesh_->cellCount(); ++cell) {
        integratedSource_[cell] = source(mesh_->cellCentroid(cell)) * mesh_->cellVolume(cell);
    }
}

std::vector<double> UnstructuredScalarTransportSolver::localTimeSteps() const {
    // **One step size per cell, not one for the mesh.** Only the steady state
    // is wanted here, and it does not depend on the path taken to it, so
    // nothing is gained by marching every cell on the clock of the smallest
    // sliver -- which is what a global step does, and these meshes carry a
    // volume ratio around 20. The transient this produces is not physical;
    // the fixed point it converges to is exactly the same one.
    //
    // The operator stays symmetric because V/dt enters only the diagonal, so
    // the base's Conjugate Gradient still applies -- which is the reason this
    // acceleration is available at all here.
    std::vector<double> outflow(mesh_->cellCount(), 0.0);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const double magnitude = std::fabs(interiorMassFlux_[i]);
        outflow[interiorFaces_[i].owner] += magnitude;
        outflow[interiorFaces_[i].neighbour] += magnitude;
    }
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        outflow[boundaryFaces_[i].cell] += std::fabs(boundaryMassFlux_[i]);
    }
    std::vector<double> steps(mesh_->cellCount(), 0.0);
    double fallback = 1e300;
    for (std::size_t cell = 0; cell < mesh_->cellCount(); ++cell) {
        if (outflow[cell] > 0.0) {
            // Safety factor, for the reason ROADMAP Fase 1 recorded: a scheme
            // sitting at CFL 1.0000 fails by going to NaN, not by warning.
            steps[cell] = 0.4 * mesh_->cellVolume(cell) / outflow[cell];
            fallback = std::min(fallback, steps[cell]);
        }
    }
    // A cell with no flux through any face has no convective limit at all;
    // give it the strictest one rather than an infinite step.
    for (double& step : steps) {
        if (step == 0.0) {
            step = fallback;
        }
    }
    return steps;
}

double UnstructuredScalarTransportSolver::stableTimeStep() const {
    const std::vector<double> steps = localTimeSteps();
    return *std::min_element(steps.begin(), steps.end());
}

double UnstructuredScalarTransportSolver::maxCellPeclet() const {
    if (diffusivity_ <= 0.0) {
        return 1e300;
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const double speed = std::fabs(interiorMassFlux_[i]) / face.areaVector.norm();
        worst = std::max(worst, speed * face.distance / diffusivity_);
    }
    return worst;
}

std::vector<double> UnstructuredScalarTransportSolver::convectiveFlux(
    const std::vector<double>& phi) const {
    // Conservative form: what leaves one cell enters its neighbour exactly,
    // because both sides use the same face value.
    const std::vector<Vector3> gradients = computeCellGradients(phi);
    std::vector<double> result(phi.size(), 0.0);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const double massFlux = interiorMassFlux_[i];
        const double flux = massFlux * faceValue(face, massFlux, phi, gradients, scheme_);
        result[face.owner] -= flux;
        result[face.neighbour] += flux;
    }
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        const BoundaryFace& face = boundaryFaces_[i];
        // A prescribed value is known exactly on the face, so it is used
        // whichever way the flow crosses -- upwinding a value that is not an
        // approximation would only throw accuracy away. An unprescribed face
        // carries the cell's own value, which is the zero-gradient outflow
        // condition.
        const double faceScalar =
            boundaryIsDirichlet_[face.meshFace] ? boundaryValue_[face.meshFace] : phi[face.cell];
        result[face.cell] -= boundaryMassFlux_[i] * faceScalar;
    }
    return result;
}

std::vector<double> UnstructuredScalarTransportSolver::applyTransportOperator(
    const std::vector<double>& x, const std::vector<double>& dt) const {
    // (V/dt) x + Gamma * L x, with L the shared Laplacian. Symmetric positive
    // definite -- the V/dt shift only adds to the diagonal -- so the base's
    // Conjugate Gradient applies unchanged. Convection is not in here: it is
    // explicit, which is what keeps this operator symmetric.
    std::vector<double> result = applyLaplacian(x);
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = mesh_->cellVolume(cell) / dt[cell] * x[cell] + diffusivity_ * result[cell];
    }
    return result;
}

std::size_t UnstructuredScalarTransportSolver::solveSteady(double tolerance, std::size_t maxSteps) {
    const std::size_t n = phi_.size();
    const std::vector<double> dt = localTimeSteps();

    // Time-independent part of the right-hand side: the volumetric source and
    // the diffusive flux from each Dirichlet face.
    std::vector<double> steadyRhs(n, 0.0);
    for (const BoundaryFace& face : boundaryFaces_) {
        if (face.entersOperator) {
            steadyRhs[face.cell] += diffusivity_ * face.coefficient * boundaryValue_[face.meshFace];
        }
    }
    if (!integratedSource_.empty()) {
        for (std::size_t cell = 0; cell < n; ++cell) {
            steadyRhs[cell] += integratedSource_[cell];
        }
    }

    std::size_t step = 0;
    for (; step < maxSteps; ++step) {
        const std::vector<double> convection = convectiveFlux(phi_);
        // The non-orthogonal part of the diffusive flux, deferred to the
        // right-hand side exactly as the steady diffusion solver does it --
        // without this the diffusion is the version Fase 2.2 measured
        // stagnating at order 0.10.
        const std::vector<double> correction = nonOrthogonalCorrection(phi_);

        std::vector<double> rhs(n);
        for (std::size_t cell = 0; cell < n; ++cell) {
            rhs[cell] = mesh_->cellVolume(cell) / dt[cell] * phi_[cell] + convection[cell] +
                        steadyRhs[cell] + diffusivity_ * correction[cell];
        }

        const std::vector<double> previous = phi_;
        conjugateGradient([this, &dt](const std::vector<double>& v) {
            return applyTransportOperator(v, dt);
        }, rhs, phi_, n, 1e-14);

        double change = 0.0;
        for (std::size_t cell = 0; cell < n; ++cell) {
            change = std::max(change, std::fabs(phi_[cell] - previous[cell]));
        }
        lastChange_ = change;
        if (change < tolerance) {
            ++step;
            break;
        }
    }
    return step;
}

} // namespace aether::solver
