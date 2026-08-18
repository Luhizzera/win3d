#include "aether/solver/UnstructuredDiffusionSolver.hpp"

#include <stdexcept>

namespace aether::solver {

using core::Vector3;
using mesh::TetrahedralMesh;

UnstructuredDiffusionSolver::UnstructuredDiffusionSolver(const TetrahedralMesh& mesh)
    : UnstructuredFvmBase(mesh) {
    phi_.assign(mesh.cellCount(), 0.0);
    boundaryIsDirichlet_.assign(mesh.faceCount(), false);
    boundaryValue_.assign(mesh.faceCount(), 0.0);
    // Registered once: the lambda reads the arrays above, so it stays correct
    // across every later rebuild.
    setBoundaryFaceValue([this](std::size_t meshFace, double& value) {
        if (!boundaryIsDirichlet_[meshFace]) {
            return false;
        }
        value = boundaryValue_[meshFace];
        return true;
    });
    rebuildCoefficients();
}

void UnstructuredDiffusionSolver::setDirichletBoundary(
    const std::function<bool(const Vector3&)>& selector, double value) {
    setDirichletBoundary(selector, [value](const Vector3&) { return value; });
}

void UnstructuredDiffusionSolver::setDirichletBoundary(
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
    rebuildCoefficients();
}

void UnstructuredDiffusionSolver::rebuildCoefficients() {
    // A Dirichlet face is an extra connection and enters the operator; an
    // insulated one carries no flux, so it contributes no coefficient at all.
    buildFaceGeometry([this](std::size_t meshFace) { return boundaryIsDirichlet_[meshFace]; });

    hasDirichletFace_ = false;
    for (const BoundaryFace& face : boundaryFaces_) {
        if (face.entersOperator) {
            hasDirichletFace_ = true;
            break;
        }
    }

    // A fixed value at the boundary is real information about the gradient in
    // the cell behind it, so it belongs in the least-squares fit.
    buildGradientStencils(/*useBoundaryValues=*/true);
}

void UnstructuredDiffusionSolver::setSourceTerm(std::function<double(const Vector3&)> source) {
    if (!source) {
        integratedSource_.clear();
        return;
    }
    integratedSource_.assign(mesh_->cellCount(), 0.0);
    for (std::size_t cell = 0; cell < mesh_->cellCount(); ++cell) {
        integratedSource_[cell] = source(mesh_->cellCentroid(cell)) * mesh_->cellVolume(cell);
    }
}

std::size_t UnstructuredDiffusionSolver::solveConjugateGradient(std::size_t maxIterations, double tolerance,
                                                                 std::size_t maxOuterSweeps) {
    if (!hasDirichletFace_) {
        throw std::runtime_error(
            "UnstructuredDiffusionSolver: needs at least one Dirichlet boundary face; with pure Neumann "
            "boundaries the operator is singular");
    }

    // Boundary source: a_b * phi_b from each face held at a fixed value.
    std::vector<double> baseRhs(phi_.size(), 0.0);
    for (const BoundaryFace& face : boundaryFaces_) {
        if (face.entersOperator) {
            baseRhs[face.cell] += face.coefficient * boundaryValue_[face.meshFace];
        }
    }
    // Volumetric source, already integrated over each cell. The operator
    // assembled here is -integral(laplacian(phi)) over the cell, so a source
    // on this side of the equation solves laplacian(phi) + S = 0.
    if (!integratedSource_.empty()) {
        for (std::size_t cell = 0; cell < phi_.size(); ++cell) {
            baseRhs[cell] += integratedSource_[cell];
        }
    }

    // Fresh each solve: this class runs the whole iteration in one call, so
    // there is nothing to carry between calls.
    double relaxation = 1.0;
    return solveDeferredCorrection(phi_, baseRhs, kNoPinnedCell, maxIterations, tolerance, maxOuterSweeps,
                                    lastOuterChange_, relaxation);
}

} // namespace aether::solver
