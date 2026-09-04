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

    // buildFaceGeometry() just rebuilt interiorFaces_/boundaryFaces_ from
    // scratch, which interiorFaceWeight_/boundaryFaceWeight_ are indexed
    // parallel to -- stale after any rebuild, not just the first one. A
    // caller is free to call setConductivity() before or after
    // setDirichletBoundary(); refreshing here whenever a conductivity is
    // already registered is what makes the order not matter.
    if (hasConductivity_) {
        updateConductivity();
    }
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

    // Boundary source: a_b * phi_b from each face held at a fixed value,
    // scaled by that face's conductivity when one was set -- the same
    // weighting applyWeightedLaplacian's diagonal uses for that face, so the
    // two sides of the equation stay consistent. weight is exactly 1.0 when
    // no conductivity was set, reproducing the plain sum below bit-for-bit.
    std::vector<double> baseRhs(phi_.size(), 0.0);
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        const BoundaryFace& face = boundaryFaces_[i];
        if (!face.entersOperator) {
            continue;
        }
        const double weight = hasConductivity_ ? boundaryFaceWeight_[i] : 1.0;
        baseRhs[face.cell] += weight * face.coefficient * boundaryValue_[face.meshFace];
    }
    // Volumetric source, already integrated over each cell. The operator
    // assembled here is -integral(laplacian(phi)) over the cell, so a source
    // on this side of the equation solves laplacian(phi) + S = 0. Not scaled
    // by conductivity: this is a heat generation rate, independent of the
    // local material's ability to conduct it away.
    if (!integratedSource_.empty()) {
        for (std::size_t cell = 0; cell < phi_.size(); ++cell) {
            baseRhs[cell] += integratedSource_[cell];
        }
    }

    // Fresh each solve: this class runs the whole iteration in one call, so
    // there is nothing to carry between calls.
    double relaxation = 1.0;
    if (hasConductivity_) {
        return solveWeightedDeferredCorrection(phi_, baseRhs, interiorFaceWeight_, boundaryFaceWeight_,
                                                kNoPinnedCell, maxIterations, tolerance, maxOuterSweeps,
                                                lastOuterChange_, relaxation);
    }
    return solveDeferredCorrection(phi_, baseRhs, kNoPinnedCell, maxIterations, tolerance, maxOuterSweeps,
                                    lastOuterChange_, relaxation);
}

void UnstructuredDiffusionSolver::setConductivity(std::function<double(const Vector3&)> conductivity) {
    if (!conductivity) {
        hasConductivity_ = false;
        return;
    }
    cellConductivity_.assign(mesh_->cellCount(), 0.0);
    for (std::size_t cell = 0; cell < mesh_->cellCount(); ++cell) {
        cellConductivity_[cell] = conductivity(mesh_->cellCentroid(cell));
    }
    hasConductivity_ = true;
    updateConductivity();
}

void UnstructuredDiffusionSolver::updateConductivity() {
    // Distance-weighted harmonic mean -- the two-point-equivalent face
    // coefficient that continuity of flux across a material interface
    // implies. ownerWeight = d_N/(d_P+d_N) (UnstructuredFvmBase.cpp), so
    // d_P = (1-ownerWeight)*distance and d_N = ownerWeight*distance; solving
    // q = k_P(T_P-T_face)/d_P = k_N(T_face-T_N)/d_N for T_face and
    // substituting back into either side gives a two-point flux
    // q = k_face(T_N-T_P)/distance with
    //
    //   k_face = 1 / ((1-ownerWeight)/k_P + ownerWeight/k_N)
    //
    // which reduces to k_face = k when k_P == k_N == k, so a uniform
    // conductivity field reproduces the unweighted operator exactly -- the
    // property that keeps setConductivity() opt-in.
    interiorFaceWeight_.resize(interiorFaces_.size());
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const double kOwner = cellConductivity_[face.owner];
        const double kNeighbour = cellConductivity_[face.neighbour];
        interiorFaceWeight_[i] =
            1.0 / ((1.0 - face.ownerWeight) / kOwner + face.ownerWeight / kNeighbour);
    }
    // A boundary face has no "other side": the prescribed value sits exactly
    // at the face, so its flux uses only the owning cell's own conductivity.
    boundaryFaceWeight_.resize(boundaryFaces_.size());
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        boundaryFaceWeight_[i] = cellConductivity_[boundaryFaces_[i].cell];
    }
}

} // namespace aether::solver
