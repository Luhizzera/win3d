#include "aether/solver/UnstructuredCavitySolver3D.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace aether::solver {

using core::Vector3;
using mesh::TetrahedralMesh;

UnstructuredCavitySolver3D::UnstructuredCavitySolver3D(
    const TetrahedralMesh& mesh, double viscosity,
    std::function<Vector3(const Vector3&)> wallVelocity,
    std::function<bool(const Vector3&)> isOutlet, double outletPressure,
    std::size_t pressureCorrectors, ConvectionScheme convection, TurbulenceModel turbulence)
    : UnstructuredFvmBase(mesh), viscosity_(viscosity), wallVelocity_(std::move(wallVelocity)),
      isOutlet_(std::move(isOutlet)), outletPressure_(outletPressure), convection_(convection),
      pressureCorrectors_(std::max<std::size_t>(pressureCorrectors, 2)), turbulence_(turbulence) {
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
    // What a boundary face's *pressure* is: the prescribed value at an
    // outlet, nothing at a wall (zero-gradient, so the face carries the
    // cell's own value). Used by the Green-Gauss fallback and by the
    // non-orthogonal correction.
    setBoundaryFaceValue([this](std::size_t meshFace, double& value) {
        if (!isOutlet_ || !isOutlet_(mesh_->face(meshFace).centroid)) {
            return false;
        }
        value = outletPressure_;
        return true;
    });
    buildBoundaryConditions();

    // Interior neighbours **and the outlet**, which is the whole content of
    // the flag. A solid wall carries a zero-gradient condition for pressure,
    // so including it would assert "the wall value equals the cell value" and
    // bias the fit -- that face carries no information about the gradient.
    // An outlet is the opposite case: its pressure is prescribed, exactly
    // like the hot and cold faces the diffusion solver has always included.
    // Leaving it out was a third way the two copies had silently diverged
    // (DIVIDA_TECNICA.md 2.3); setBoundaryFaceValue() above already answers
    // "does this face have a value", and only an outlet says yes.
    buildGradientStencils(/*useBoundaryValues=*/true);
    boundaryFlux_.assign(boundaryFaces_.size(), 0.0);

    eddyViscosity_.assign(mesh.cellCount(), 0.0);
    if (turbulence_ != TurbulenceModel::None) {
        buildWallDistances();
    }
    updateEffectiveViscosity();
}

void UnstructuredCavitySolver3D::buildWallDistances() {
    wallDistance_.assign(mesh_->cellCount(), 0.0);
    for (std::size_t cell = 0; cell < mesh_->cellCount(); ++cell) {
        const Vector3& centroid = mesh_->cellCentroid(cell);
        double best = std::numeric_limits<double>::max();
        for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
            if (boundaryConditions_[i].isOutlet) {
                continue; // an outlet is not a wall: nothing damps turbulence there
            }
            best = std::min(best, (centroid - boundaryFaces_[i].centroid).norm());
        }
        // A domain with no solid wall at all leaves `best` at its sentinel;
        // the cap below is what bounds the mixing length in that case, and
        // an unbounded wall distance would otherwise make l_m meaningless.
        wallDistance_[cell] = best == std::numeric_limits<double>::max() ? 0.0 : best;
    }

    // The mixing length's far-field cap, from the domain's own smallest
    // extent -- the same 0.09 * (half the smallest side) that every
    // mixing-length closure in this project uses (the Escudier asymptotic
    // factor). Computed from the mesh bounding box because an unstructured
    // domain has no Lx/Ly/Lz to be told.
    Vector3 lo = mesh_->vertex(0);
    Vector3 hi = lo;
    for (std::size_t v = 1; v < mesh_->vertexCount(); ++v) {
        const Vector3& q = mesh_->vertex(v);
        lo = Vector3(std::min(lo.x, q.x), std::min(lo.y, q.y), std::min(lo.z, q.z));
        hi = Vector3(std::max(hi.x, q.x), std::max(hi.y, q.y), std::max(hi.z, q.z));
    }
    const double smallestExtent = std::min({hi.x - lo.x, hi.y - lo.y, hi.z - lo.z});
    mixingLengthCap_ = 0.09 * smallestExtent * 0.5;
}

void UnstructuredCavitySolver3D::updateEddyViscosity(
    const std::vector<std::vector<Vector3>>& velocityGradient) {
    if (turbulence_ == TurbulenceModel::None) {
        return;
    }
    constexpr double kVonKarman = 0.41;
    for (std::size_t cell = 0; cell < mesh_->cellCount(); ++cell) {
        const Vector3& gu = velocityGradient[0][cell];
        const Vector3& gv = velocityGradient[1][cell];
        const Vector3& gw = velocityGradient[2][cell];
        // S_ij = (du_i/dx_j + du_j/dx_i)/2, and |S| = sqrt(2 S_ij S_ij).
        // Written out rather than looped because the six independent
        // components each name a different pair of gradients, and an
        // indexed form obscures which is which.
        const double sxx = gu.x;
        const double syy = gv.y;
        const double szz = gw.z;
        const double sxy = 0.5 * (gu.y + gv.x);
        const double sxz = 0.5 * (gu.z + gw.x);
        const double syz = 0.5 * (gv.z + gw.y);
        const double strain = std::sqrt(2.0 * (sxx * sxx + syy * syy + szz * szz +
                                                2.0 * (sxy * sxy + sxz * sxz + syz * syz)));
        const double mixingLength = std::min(kVonKarman * wallDistance_[cell], mixingLengthCap_);
        eddyViscosity_[cell] = mixingLength * mixingLength * strain;
    }
    updateEffectiveViscosity();
}

void UnstructuredCavitySolver3D::updateEffectiveViscosity() {
    const std::size_t n = mesh_->cellCount();
    viscousDiagonal_.assign(n, 0.0);

    if (turbulence_ == TurbulenceModel::None) {
        // Exactly the expression the laminar path always used, so its
        // results stay bit-identical rather than merely close -- see this
        // function's own declaration for why that is worth a branch.
        interiorFaceViscosity_.assign(interiorFaces_.size(), viscosity_);
        boundaryFaceViscosity_.assign(boundaryFaces_.size(), viscosity_);
        for (std::size_t cell = 0; cell < n; ++cell) {
            viscousDiagonal_[cell] = viscosity_ * interiorDiagonal_[cell];
        }
        return;
    }

    interiorFaceViscosity_.resize(interiorFaces_.size());
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const double nu =
            viscosity_ + 0.5 * (eddyViscosity_[face.owner] + eddyViscosity_[face.neighbour]);
        interiorFaceViscosity_[i] = nu;
        viscousDiagonal_[face.owner] += nu * face.coefficient;
        viscousDiagonal_[face.neighbour] += nu * face.coefficient;
    }
    // **Molecular viscosity only on a solid wall face**, not the adjacent
    // cell's turbulent value. The mixing length -- and therefore nu_t -- is
    // exactly zero at a wall by definition, and this project already paid
    // for getting that wrong once: MixingLengthChannelFlowSolver1D mirrored
    // the neighbouring cell's nu_t across the wall and its friction
    // velocity came out 38% off the exact momentum balance.
    boundaryFaceViscosity_.assign(boundaryFaces_.size(), viscosity_);
}

void UnstructuredCavitySolver3D::loadState(std::vector<Vector3> velocity, std::vector<double> pressure,
                                            double time) {
    if (velocity.size() != velocity_.size() || pressure.size() != pressure_.size()) {
        throw std::invalid_argument(
            "UnstructuredCavitySolver3D::loadState: field size does not match the mesh's cell count");
    }
    velocity_ = std::move(velocity);
    pressure_ = std::move(pressure);
    time_ = time;
    // The recorded boundary flux belongs to the projection that produced the
    // old field, not to this one; leaving it would make the first
    // netBoundaryFlux() after a load report the previous run's outlets.
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
                                                                 double dt, bool correctionField) const {
    // Rhie-Chow style: interpolate the velocity to the face, then replace the
    // pressure-gradient part of that interpolation with the compact face
    // difference. Without this the face flux is blind to a cell-to-cell
    // pressure oscillation, which is the collocated-grid failure mode this
    // arrangement is otherwise exposed to.
    const std::vector<Vector3> pressureGradients =
        correctionField ? computeCellGradients(pressure, /*clampFallback=*/false,
                                                /*homogeneousBoundaryValues=*/true)
                        : computeCellGradients(pressure);
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

std::vector<double> UnstructuredCavitySolver3D::applyHelmholtzOperator(
    const std::vector<double>& x, double dt, const std::vector<double>& convectionOutflow,
    bool viscous) const {
    std::vector<double> result(x.size());
    for (std::size_t cell = 0; cell < x.size(); ++cell) {
        result[cell] = (mesh_->cellVolume(cell) / dt +
                        (viscous ? viscousDiagonal_[cell] : 0.0) + convectionOutflow[cell]) *
                       x[cell];
    }
    if (!viscous) {
        return result;
    }
    // A wall face stiffens its cell's equation even though the wall value
    // itself is known and belongs on the right-hand side.
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        if (boundaryConditions_[i].isOutlet) {
            continue; // zero-gradient: no viscous flux through an outlet
        }
        result[boundaryFaces_[i].cell] +=
            boundaryFaceViscosity_[i] * boundaryFaces_[i].coefficient * x[boundaryFaces_[i].cell];
    }
    // Still symmetric with a variable coefficient: both sides of a face use
    // the *same* face viscosity, so the off-diagonal pair stays equal and
    // the Conjugate Gradient above applies unchanged.
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const double nu = interiorFaceViscosity_[i];
        result[face.owner] -= nu * face.coefficient * x[face.neighbour];
        result[face.neighbour] -= nu * face.coefficient * x[face.owner];
    }
    return result;
}

std::vector<double> UnstructuredCavitySolver3D::solveHelmholtz(
    const std::vector<double>& rhs, double dt, const std::vector<double>& convectionOutflow,
    bool viscous) const {
    // Plain CG. This operator is far more diagonally dominant than the
    // pressure Poisson one -- the V/dt shift only adds to the diagonal -- so
    // it converges in very few iterations, which is why making diffusion
    // implicit costs much less than the step limit it removes.
    std::vector<double> x(rhs.size(), 0.0);
    conjugateGradient(
        [this, dt, &convectionOutflow, viscous](const std::vector<double>& v) {
            return applyHelmholtzOperator(v, dt, convectionOutflow, viscous);
        },
        rhs, x, rhs.size(), 1e-12);
    return x;
}

std::vector<double> UnstructuredCavitySolver3D::divergenceOfFlux(const std::vector<Vector3>& velocity,
                                                                   const std::vector<double>& pressure,
                                                                   double fluxDt, double outerDt,
                                                                   bool correctionField) const {
    const std::size_t n = velocity.size();
    const std::vector<double> fluxes = faceMassFluxes(velocity, pressure, fluxDt, correctionField);
    std::vector<double> rhs(n, 0.0);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        rhs[interiorFaces_[i].owner] -= fluxes[i] / outerDt;
        rhs[interiorFaces_[i].neighbour] += fluxes[i] / outerDt;
    }
    // Boundary contributions: a prescribed-velocity face (wall or inlet)
    // injects its known flux, and an outlet both carries flux and pins the
    // pressure level through its Dirichlet coefficient. Neither carries a
    // Rhie-Chow term at all -- fluxDt only ever affects interior faces --
    // so this half is exactly the same whatever fluxDt is.
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        const BoundaryFace& face = boundaryFaces_[i];
        const BoundaryCondition& condition = boundaryConditions_[i];
        const double flux = condition.isOutlet ? velocity[face.cell].dot(face.areaVector)
                                                : condition.wallVelocity.dot(face.areaVector);
        rhs[face.cell] -= flux / outerDt;
        if (condition.isOutlet) {
            rhs[face.cell] += face.coefficient * outletPressure_;
        }
    }
    return rhs;
}

void UnstructuredCavitySolver3D::projectToDivergenceFree(std::vector<Vector3>& velocityStar, double dt) {
    const std::size_t n = velocityStar.size();
    // Cell 0 is pinned ONLY when every boundary is a wall: with an outlet the
    // prescribed pressure already fixes the level, and pinning as well would
    // over-constrain the system with two incompatible references.
    const std::size_t pinnedCell = hasOutlet_ ? kNoPinnedCell : 0;

    // **First corrector: unchanged from before this item, bit-identical.**
    // The same Laplacian the diffusion solver uses, correction included --
    // this operator used to be a second copy that assembled only the
    // implicit part, the version Fase 2.2 measured stagnating at observed
    // order 0.10 on exactly these meshes, in the most important equation
    // this solver has (DIVIDA_TECNICA.md 1.2). fluxDt=0 here is deliberate:
    // the predictor carries no pressure correction yet, so there is no
    // Rhie-Chow term to be consistent with. Started from the previous step's
    // pressure, which is why it converges in a handful of iterations.
    std::vector<double> rhs = divergenceOfFlux(velocityStar, pressure_, 0.0, dt);
    solveDeferredCorrection(pressure_, rhs, pinnedCell, n, 1e-10, pressureCorrectors_,
                            lastPressureChange_, pressureRelaxation_);

    std::vector<Vector3> pressureGradients = computeCellGradients(pressure_);
    for (std::size_t cell = 0; cell < n; ++cell) {
        velocity_[cell] = velocityStar[cell] - pressureGradients[cell] * dt;
    }

    // Record the boundary flux **the first corrector enforced**, term for
    // term, *before* the coupling correction below touches pressure_ --
    // this exact formula's mass-conservation guarantee (measured to close
    // the channel to 6e-14) is an algebraic identity of *this* Laplacian
    // system specifically (its RHS's own boundary terms telescope against
    // this face flux), and does not carry over to whatever pressure_ becomes
    // after a second, different operator perturbs it. See the coupling
    // correction below for how its own contribution is added back in a way
    // that *is* consistent with the operator that produced it.
    //
    // Both terms of the operator have to appear here, and that was measured:
    // when the non-orthogonal correction was added to the Poisson equation
    // and this flux still carried only the compact part a_b (p_out - p_P),
    // the imbalance went from 6e-14 straight back to 2.5e-03. The projection
    // had simply started enforcing a flux with one more term in it. Same
    // lesson as ROADMAP Fase 1 and as item 1.1, for the third time:
    // **measure the operator you actually solve.**
    //
    // Written as (grad p)_f . A rather than as its two pieces because
    // that identity is the reason the split is legitimate at all:
    // A = a_b d + A_nonorth and (grad p)_f . d = p_out - p_P, so
    // (grad p)_f . A  ==  a_b (p_out - p_P) + (grad p)_f . A_nonorth.
    boundaryFlux_.assign(boundaryFaces_.size(), 0.0);
    for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
        const BoundaryFace& face = boundaryFaces_[i];
        if (!boundaryConditions_[i].isOutlet) {
            continue;
        }
        const Vector3& cellGradient = pressureGradients[face.cell];
        const double compactNormalDerivative =
            (outletPressure_ - pressure_[face.cell]) / face.distance;
        const Vector3 faceGradient =
            cellGradient + face.unitD * (compactNormalDerivative - cellGradient.dot(face.unitD));
        boundaryFlux_[i] = velocityStar[face.cell].dot(face.areaVector) -
                            dt * faceGradient.dot(face.areaVector);
    }

    // **Coupling correction: DIVIDA_TECNICA.md 4.3 -- still gated to closed
    // domains (no outlet), now for a narrower and better-understood
    // reason.** The sixth attempt fixed a genuine defect in this operator:
    // `computeCellGradients` was feeding outlet-adjacent cells the
    // *prescribed* `outletPressure_` as the boundary neighbour's value even
    // when the field being differentiated is a pressure *correction*, whose
    // value at a Dirichlet boundary is necessarily zero. That made the
    // operator affine rather than linear there, which is the one thing a
    // Krylov method cannot tolerate. The fix (`homogeneousBoundaryValues`)
    // is kept because it is correct independently of what it repairs, and
    // it is a no-op on a closed domain, where no boundary face carries a
    // prescribed pressure at all.
    //
    // **It improved the outlet case and did not close it**, measured on the
    // channel: the per-step residual after correction fell from 0.047 to
    // 0.033 (~30%). Enabling the correction on an outlet domain in that
    // state is still a regression rather than a partial win -- the
    // channel's mass imbalance goes from 6.5e-14 (correction skipped) to
    // 1.2e-4 (correction applied) -- so the gate stays.
    //
    // **Why it does not close is now known, and it is not what an earlier
    // version of this comment said.** That version claimed GMRES suffered
    // a numerical breakdown. Instrumenting the two breakdown paths directly
    // showed neither ever fires: GMRES *converges*, to its own 1e-8
    // relative tolerance, in 64 iterations. So the operator it solves is
    // simply not the exact linearization of the quantity that is measured
    // afterwards -- and on a boundary face those two provably disagree:
    //
    //   divergenceOfFlux (what GMRES drives to zero) uses, at an outlet,
    //     velocity[cell] . areaVector
    //   maxFaceDivergence (what actually reports mass conservation) uses
    //     boundaryFlux_[i], which is the gradient-blended formula recorded
    //     further down this function.
    //
    // That is the same class of error the fourth attempt found between
    // fluxDt=0 and fluxDt=dt, one level out: **measure the operator you
    // actually solve**, now for the third time in this item, on boundary
    // faces instead of interior ones. Reconciling the two outlet flux
    // definitions is the next concrete step, and it is a narrower question
    // than any previous framing of this gap.
    if (pinnedCell != kNoPinnedCell) {
    // The first corrector drives the fluxDt=0 residual above to ~1e-10 -- that
    // is its own convergence criterion -- but maxFaceDivergence(), the
    // metric that actually decides whether this mesh is stable, uses
    // faceMassFluxes(velocity_, pressure_, lastDt_): the **real** step dt in
    // the Rhie-Chow term. The fourth attempt built a GMRES correction
    // against the fluxDt=0 residual, the wrong target -- it solved that
    // residual beautifully (measured: ~0.08 to ~1e-9 per step) while the
    // real divergence barely moved, because the first corrector had already
    // solved that exact residual to begin with. This attempt targets the
    // real one directly: fluxDt=dt below, matching lastDt_ exactly (see
    // step(), which sets lastDt_ = dt right after this function returns).
    //
    // `correctionOperator` is the *linear* part of
    // divergenceOfFlux(·, ·, dt, dt) as a function of a pressure correction
    // x applied as an additional velocity correction -dt*grad(x): both the
    // interpolated-velocity term and the Rhie-Chow correction term inside
    // faceMassFluxes are linear in their own argument (velocity, pressure
    // respectively) at fixed fluxDt, so passing (-dt*grad(x), x) together
    // captures the *whole* linear response, not just the velocity half --
    // the fourth attempt's operator only ever varied the velocity argument
    // and left the pressure argument at pressure_, which silently dropped
    // the fluxDt-dependent term's own contribution to the linearization.
    // Boundary terms are affine, not linear (a wall's prescribed flux and an
    // outlet's `coefficient * outletPressure_` do not depend on the
    // argument at all), so that constant is subtracted off by evaluating the
    // same function at (zero velocity, zero pressure), which is exactly
    // divergenceOfFlux()'s own documented affine part.
    //
    // computeCellGradients(x, false): the *unclamped* gradient, for the same
    // reason the fourth attempt's second bug required it -- the clamped
    // version is not linear on deficient-stencil cells, and a Krylov method
    // handed an operator that lies about linearity converges confidently to
    // the wrong answer. The velocity update after solving uses the same
    // unclamped operator for this delta, incrementally on top of the first
    // corrector's velocity_, rather than a fresh clamped gradient of the new
    // total pressure -- keeping the applied correction identical to what
    // this operator promised GMRES.
    const std::vector<Vector3> zeroVelocity(n, Vector3{});
    const std::vector<double> zeroPressure(n, 0.0);
    const std::vector<double> wallOnlyDivergence =
        divergenceOfFlux(zeroVelocity, zeroPressure, dt, dt, /*correctionField=*/true);
    const auto correctionOperator = [&](const std::vector<double>& x) {
        const std::vector<Vector3> gradient =
            computeCellGradients(x, /*clampFallback=*/false, /*homogeneousBoundaryValues=*/true);
        std::vector<Vector3> correctionVelocity(n);
        for (std::size_t cell = 0; cell < n; ++cell) {
            correctionVelocity[cell] = gradient[cell] * (-dt);
        }
        std::vector<double> divergence =
            divergenceOfFlux(correctionVelocity, x, dt, dt, /*correctionField=*/true);
        for (std::size_t cell = 0; cell < n; ++cell) {
            divergence[cell] -= wallOnlyDivergence[cell];
        }
        if (pinnedCell != kNoPinnedCell) {
            divergence[pinnedCell] = x[pinnedCell]; // identity row, matching the CG operator's pinning
        }
        return divergence;
    };

    std::vector<double> residual = divergenceOfFlux(velocity_, pressure_, dt, dt);
    if (pinnedCell != kNoPinnedCell) {
        residual[pinnedCell] = 0.0;
    }
    double residualNorm = 0.0;
    for (std::size_t cell = 0; cell < n; ++cell) {
        residualNorm = std::max(residualNorm, std::fabs(residual[cell]));
    }
    lastCouplingResidualBefore_ = residualNorm;

    // Relative tolerance, learned from DIVIDA_TECNICA.md 5.4: an absolute one
    // is either too loose to matter or too tight to ever fire, depending on
    // the field's own scale. Skipped entirely below the floor rather than
    // running a GMRES cycle whose own first residual check would immediately
    // exit anyway -- cheap on the meshes that do not need it.
    const double couplingTolerance = 1e-8;
    if (residualNorm > 1e-12) {
        // Solved for -residual, not residual: correctionOperator(x) is the
        // divergence *added* by the correction, so cancelling the existing
        // residual means correctionOperator(dp) must equal its negative.
        std::vector<double> negatedResidual(n);
        for (std::size_t cell = 0; cell < n; ++cell) {
            negatedResidual[cell] = -residual[cell];
        }
        std::vector<double> deltaPressure(n, 0.0);
        // A restart this small was measured not converging within budget on
        // meshes of a few hundred cells -- GMRES(30) throws away the Krylov
        // subspace every 30 iterations, so a system that needs more than 30
        // directions to resolve well restarts from scratch repeatedly and
        // stalls well short of the requested tolerance. Unrestarted GMRES
        // converges in at most n iterations in exact arithmetic, so the
        // restart is capped at n rather than a small fixed constant.
        const std::size_t restart = std::min<std::size_t>(200, n);
        const std::size_t maxIterations = std::min<std::size_t>(2000, 4 * n);
        lastCouplingIterations_ = gmres(correctionOperator, negatedResidual, deltaPressure, restart,
                                        maxIterations, couplingTolerance);

        for (std::size_t cell = 0; cell < n; ++cell) {
            pressure_[cell] += deltaPressure[cell];
        }
        // deltaVelocity is exactly the correction correctionOperator's own
        // outlet term is linear in (see divergenceOfFlux's outlet handling:
        // plain velocity.dot(area), no gradient blend) -- adding its outlet
        // contribution to boundaryFlux_ here, on top of the first
        // corrector's already-recorded flux above, is what keeps the
        // recorded flux consistent with *whichever* operator actually moved
        // the pressure at that face, instead of re-deriving a single
        // snapshot formula that only one of the two operators satisfies.
        const std::vector<Vector3> deltaGradient = computeCellGradients(
            deltaPressure, /*clampFallback=*/false, /*homogeneousBoundaryValues=*/true);
        std::vector<Vector3> deltaVelocity(n);
        for (std::size_t cell = 0; cell < n; ++cell) {
            deltaVelocity[cell] = deltaGradient[cell] * (-dt);
            velocity_[cell] += deltaVelocity[cell];
        }
        for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
            if (!boundaryConditions_[i].isOutlet) {
                continue;
            }
            boundaryFlux_[i] += deltaVelocity[boundaryFaces_[i].cell].dot(boundaryFaces_[i].areaVector);
        }

        std::vector<double> finalResidual = divergenceOfFlux(velocity_, pressure_, dt, dt);
        if (pinnedCell != kNoPinnedCell) {
            finalResidual[pinnedCell] = 0.0;
        }
        double finalNorm = 0.0;
        for (std::size_t cell = 0; cell < n; ++cell) {
            finalNorm = std::max(finalNorm, std::fabs(finalResidual[cell]));
        }
        lastCouplingResidualAfter_ = finalNorm;
    } else {
        lastCouplingIterations_ = 0;
        lastCouplingResidualAfter_ = residualNorm;
    }
    } else {
        // Outlet mesh: correction skipped, see above. Zeroed rather than
        // left stale, so a caller cannot mistake a previous step's numbers
        // for this one's.
        lastCouplingIterations_ = 0;
        lastCouplingResidualBefore_ = 0.0;
        lastCouplingResidualAfter_ = 0.0;
    }
}

void UnstructuredCavitySolver3D::step(double dt) { stepWith(dt, StepParts{}); }

void UnstructuredCavitySolver3D::stepWith(double dt, StepParts parts) {
    const std::size_t n = velocity_.size();
    std::vector<Vector3> flux(n, Vector3{});

    // --- Convection, upwinded on the face mass flux (see the class comment
    // for why not central differencing), in conservative form so that what
    // leaves one cell enters its neighbour exactly.
    const std::vector<double> massFluxes = faceMassFluxes(velocity_, pressure_, lastDt_);
    // The face reconstruction is per component, so the three components need
    // their own gradients -- the limiter is built from the upwind cell's
    // gradient of the quantity being convected, and momentum has three of
    // them. Three extra gradient passes per step, against several matrix-free
    // CG solves in the same step.
    std::vector<std::vector<double>> velocityComponent(3, std::vector<double>(n));
    for (std::size_t cell = 0; cell < n; ++cell) {
        velocityComponent[0][cell] = velocity_[cell].x;
        velocityComponent[1][cell] = velocity_[cell].y;
        velocityComponent[2][cell] = velocity_[cell].z;
    }
    std::vector<std::vector<Vector3>> velocityGradient(3);
    for (int c = 0; c < 3; ++c) {
        velocityGradient[c] = computeCellGradients(velocityComponent[c]);
    }
    // The turbulence closure rides along on the gradients the convection
    // limiter already needed, so it adds no gradient pass of its own. A
    // no-op when the closure is None, which is why it is unconditional.
    updateEddyViscosity(velocityGradient);
    // **The outflow part of convection goes implicit**, which is what removes
    // the convective step limit (DIVIDA_TECNICA.md 4.2). Accumulated here as a
    // per-cell coefficient: the mass flux leaving through every face where
    // this cell is the upwind one.
    //
    // Why this particular split, and why it keeps the Conjugate Gradient this
    // project already has: a fully implicit upwind operator is
    // **non-symmetric** -- the two neighbour coefficients differ by flow
    // direction -- and would need BiCGSTAB or GMRES. Moving only the outflow
    // makes the addition purely diagonal, so the operator stays symmetric
    // positive definite and nothing about the linear solve changes.
    //
    // And it is unconditionally bounded, which is the point. Writing the
    // first-order part of the update out,
    //
    //   u* = [ (V/dt) u^n + sum_in F u_N^n ] / [ V/dt + sum_out F ]
    //
    // and using sum_in F = sum_out F for a divergence-free face flux, u* is a
    // **convex combination** of u^n and the neighbours' u^n at any dt. No step
    // size can make it overshoot.
    std::vector<double> convectionOutflow(n, 0.0);
    for (std::size_t i = 0; i < interiorFaces_.size(); ++i) {
        const InteriorFace& face = interiorFaces_[i];
        const double massFlux = massFluxes[i];
        const Vector3 faceVelocity{
            faceValue(face, massFlux, velocityComponent[0], velocityGradient[0], convection_),
            faceValue(face, massFlux, velocityComponent[1], velocityGradient[1], convection_),
            faceValue(face, massFlux, velocityComponent[2], velocityGradient[2], convection_)};
        flux[face.owner] -= faceVelocity * massFlux;
        flux[face.neighbour] += faceVelocity * massFlux;
        convectionOutflow[massFlux >= 0.0 ? face.owner : face.neighbour] += std::fabs(massFlux);
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
        if (massFlux >= 0.0) {
            // Leaving through this face carries the cell value away, the same
            // implicit outflow term an interior face contributes. Inflow
            // carries a *prescribed* value and stays a source.
            convectionOutflow[face.cell] += massFlux;
        }
    }

    if (!parts.convection) {
        std::fill(flux.begin(), flux.end(), Vector3{});
        std::fill(convectionOutflow.begin(), convectionOutflow.end(), 0.0);
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
        // Adding the outflow term back cancels, exactly, the part of `flux`
        // that carried this cell own value away: on a face where this cell is
        // upwind, the upwind component of the face value *is* u^n here. What
        // stays explicit is only the limiter correction to it -- the same
        // deferred-correction structure the non-orthogonal term uses. The
        // cancelled part reappears on the operator diagonal.
        component[0][cell] += convectionOutflow[cell] * velocity_[cell].x;
        component[1][cell] += convectionOutflow[cell] * velocity_[cell].y;
        component[2][cell] += convectionOutflow[cell] * velocity_[cell].z;
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

    if (!parts.viscous) {
        // Undo the wall viscous source added just above; the operator's own
        // viscous terms are switched off inside solveHelmholtz by the same
        // flag, so the two have to be dropped together or the equation is
        // inconsistent.
        for (std::size_t i = 0; i < boundaryFaces_.size(); ++i) {
            const BoundaryCondition& condition = boundaryConditions_[i];
            if (condition.isOutlet) {
                continue;
            }
            const std::size_t cell = boundaryFaces_[i].cell;
            const double coefficient = viscosity_ * boundaryFaces_[i].coefficient;
            component[0][cell] -= coefficient * condition.wallVelocity.x;
            component[1][cell] -= coefficient * condition.wallVelocity.y;
            component[2][cell] -= coefficient * condition.wallVelocity.z;
        }
    }

    std::vector<Vector3> velocityStar(n);
    for (int c = 0; c < 3; ++c) {
        const std::vector<double> solved =
            solveHelmholtz(component[c], dt, convectionOutflow, parts.viscous);
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

    if (parts.projection) {
        projectToDivergenceFree(velocityStar, dt);
    } else {
        velocity_ = velocityStar;
    }
    lastDt_ = dt;
    time_ += dt;

    // **Refuse rather than propagate.** On a sufficiently distorted mesh this
    // scheme is linearly unstable and the field runs away to inf and then to
    // NaN (DIVIDA_TECNICA.md 4.3). Every number this class returns afterwards
    // is NaN, including the diagnostics a caller would check -- so without
    // this the failure is silent in the worst way: a whole run of plausible
    // machinery producing a field that is not a field.
    //
    // The check is here rather than in a mesh-quality gate at construction
    // because **there is no honest a-priori criterion to gate on**, and that
    // was measured, not assumed: a mesh with non-orthogonality 2.24 runs
    // fine, while one at 2.07 diverges. Whatever predicts this failure, the
    // quantities this class currently computes do not.
    //
    // Costs one pass over the cells per step, against several matrix-free CG
    // solves in the same step.
    for (std::size_t cell = 0; cell < n; ++cell) {
        if (std::isfinite(pressure_[cell]) && std::isfinite(velocity_[cell].x) &&
            std::isfinite(velocity_[cell].y) && std::isfinite(velocity_[cell].z)) {
            continue;
        }
        throw std::runtime_error(
            "UnstructuredCavitySolver3D: the solution stopped being finite. This is a linear "
            "instability of the scheme on this mesh, not a time-step problem: it was measured to "
            "grow by a fixed factor per step, unchanged by a 10x smaller dt or a 1e-6 smaller "
            "velocity scale, and it appears on meshes whose non-orthogonality is no worse than "
            "meshes that run fine. See DIVIDA_TECNICA.md 4.3.");
    }
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
