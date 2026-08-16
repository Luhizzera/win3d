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

std::size_t UnstructuredDiffusionSolver::solveConjugateGradient(std::size_t maxIterations, double tolerance,
                                                                 std::size_t maxOuterSweeps) {
    if (!hasDirichletFace_) {
        throw std::runtime_error(
            "UnstructuredDiffusionSolver: needs at least one Dirichlet boundary face; with pure Neumann "
            "boundaries the operator is singular");
    }

    // The only source term: a_b * phi_b from each face held at a fixed value.
    std::vector<double> baseRhs(phi_.size(), 0.0);
    for (const BoundaryFace& face : boundaryFaces_) {
        if (face.entersOperator) {
            baseRhs[face.cell] += face.coefficient * boundaryValue_[face.meshFace];
        }
    }

    return solveDeferredCorrection(phi_, baseRhs, kNoPinnedCell, maxIterations, tolerance, maxOuterSweeps,
                                    lastOuterChange_);
}

} // namespace aether::solver
