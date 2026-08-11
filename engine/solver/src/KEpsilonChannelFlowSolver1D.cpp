#include "aether/solver/KEpsilonChannelFlowSolver1D.hpp"

#include "aether/solver/MixingLengthChannelFlowSolver1D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

namespace {
constexpr double kKarman = 0.41;
constexpr double kCmu = 0.09;
constexpr double kCeps1 = 1.44;
constexpr double kCeps2 = 1.92;
constexpr double kSigmaK = 1.0;
constexpr double kSigmaEps = 1.3;
constexpr double kLogLawB = 5.0;  // standard log-law additive constant
constexpr double kFloor = 1e-12;  // guards against division by ~0 during transient Picard iterates
} // namespace

KEpsilonChannelFlowSolver1D::KEpsilonChannelFlowSolver1D(std::size_t ny, double height,
                                                           double kinematicViscosity, double source)
    : ny_(ny), height_(height), h_(height / static_cast<double>(ny)), nu_(kinematicViscosity),
      source_(source), u_(ny, 0.0), k_(ny, 0.0), epsilon_(ny, 0.0), nut_(ny, 0.0), production_(ny, 0.0) {
    const double uTau = frictionVelocity();
    const double kWall = uTau * uTau / std::sqrt(kCmu);
    for (std::size_t j = 0; j < ny_; ++j) {
        k_[j] = kWall;
        epsilon_[j] = uTau * uTau * uTau / (kKarman * std::max(wallDistance(j), h_ * 0.5));
    }

    // Warm-start the velocity field -- see the class comment for why
    // starting from rest deadlocks the k-epsilon coupling.
    MixingLengthChannelFlowSolver1D primer(ny, height, kinematicViscosity, source);
    primer.solve();
    for (std::size_t j = 0; j < ny_; ++j) {
        u_[j] = primer.u(j);
    }
}

double KEpsilonChannelFlowSolver1D::cellCenterY(std::size_t j) const {
    return (static_cast<double>(j) + 0.5) * h_;
}

double KEpsilonChannelFlowSolver1D::wallDistance(std::size_t j) const {
    const double y = cellCenterY(j);
    return std::min(y, height_ - y);
}

double KEpsilonChannelFlowSolver1D::frictionVelocity() const { return std::sqrt(source_ * height_ / 2.0); }

double KEpsilonChannelFlowSolver1D::velocityGradientAt(std::size_t j) const {
    return (u_[j + 1] - u_[j - 1]) / (2.0 * h_);
}

void KEpsilonChannelFlowSolver1D::updateEddyViscosityAndProduction() {
    // Under-relaxed: nu_t = Cmu*k^2/epsilon amplifies noise in k and
    // epsilon (a small change in epsilon, especially as it trends toward
    // small values, swings nu_t sharply), and nu_t in turn controls both
    // production and the momentum equation's own diffusivity -- without
    // relaxation this feedback loop oscillates/collapses instead of
    // settling. Blending with the previous nu_t is standard practice in
    // production k-epsilon solvers (e.g. OpenFOAM's field relaxation
    // factors), not a new modeling assumption.
    constexpr double kRelaxation = 0.3;
    for (std::size_t j = 0; j < ny_; ++j) {
        const double kSafe = std::max(k_[j], kFloor);
        const double epsSafe = std::max(epsilon_[j], kFloor);
        const double nutComputed = kCmu * kSafe * kSafe / epsSafe;
        nut_[j] = kRelaxation * nutComputed + (1.0 - kRelaxation) * nut_[j];
    }
    for (std::size_t j = 1; j + 1 < ny_; ++j) {
        const double gradient = velocityGradientAt(j);
        production_[j] = nut_[j] * gradient * gradient;
    }
}

std::size_t KEpsilonChannelFlowSolver1D::solve(std::size_t maxOuterIterations, double tolerance) {
    const double uTau = frictionVelocity();
    const double kWall = uTau * uTau / std::sqrt(kCmu);
    const std::size_t last = ny_ - 1;
    const double epsWallBottom = uTau * uTau * uTau / (kKarman * wallDistance(0));
    const double epsWallTop = uTau * uTau * uTau / (kKarman * wallDistance(last));

    k_[0] = kWall;
    k_[last] = kWall;
    epsilon_[0] = epsWallBottom;
    epsilon_[last] = epsWallTop;

    // The momentum equation's wall-adjacent cells are *fixed* directly to
    // the log-law value (not solved via a wall-shear flux condition): with
    // a flux/Neumann condition at *both* walls, the momentum equation has
    // no Dirichlet anchor at all, leaving u determined only up to an
    // arbitrary additive constant (the same null-space issue the pressure
    // solves in TaylorGreenVortexSolver2D/LidDrivenCavitySolver2D have,
    // fixed there by pinning one cell) -- caught here because the solved
    // profile came out uniformly, unphysically negative despite otherwise
    // sensible k/epsilon/nu_t. Using the log law to fix these two cells
    // pulls in the empirical additive constant B (~5.0, the one part of
    // this project's turbulence validation that *is* a recalled literature
    // value) purely to anchor the absolute level; the log-law *slope*
    // check in the test does not depend on having B exactly right.
    const double yPlusBottom = wallDistance(0) * uTau / nu_;
    const double yPlusTop = wallDistance(last) * uTau / nu_;
    u_[0] = uTau * ((1.0 / kKarman) * std::log(yPlusBottom) + kLogLawB);
    u_[last] = uTau * ((1.0 / kKarman) * std::log(yPlusTop) + kLogLawB);

    std::size_t outer = 0;
    for (; outer < maxOuterIterations; ++outer) {
        updateEddyViscosityAndProduction();

        double maxDelta = 0.0;

        // Momentum: interior cells only; wall cells (0, last) stay fixed
        // at the log-law value set above (the mesh does not resolve down
        // to the wall itself).
        for (int inner = 0; inner < 20; ++inner) {
            for (std::size_t j = 1; j + 1 < ny_; ++j) {
                const double gammaHere = nu_ + nut_[j];
                const double gammaUp = 0.5 * (gammaHere + (nu_ + nut_[j + 1]));
                const double gammaDown = 0.5 * (gammaHere + (nu_ + nut_[j - 1]));
                const double newU = (gammaUp * u_[j + 1] + gammaDown * u_[j - 1] + source_ * h_ * h_) /
                                     (gammaUp + gammaDown);
                maxDelta = std::max(maxDelta, std::fabs(newU - u_[j]));
                u_[j] = newU;
            }
        }

        // k: interior cells solved; wall cells fixed at the equilibrium
        // wall-function value.
        for (int inner = 0; inner < 20; ++inner) {
            for (std::size_t j = 1; j + 1 < ny_; ++j) {
                const double gammaHere = nu_ + nut_[j] / kSigmaK;
                const double gammaUp = 0.5 * (gammaHere + (nu_ + nut_[j + 1] / kSigmaK));
                const double gammaDown = 0.5 * (gammaHere + (nu_ + nut_[j - 1] / kSigmaK));
                const double sourceTerm = (production_[j] - epsilon_[j]) * h_ * h_;
                const double newK =
                    (gammaUp * k_[j + 1] + gammaDown * k_[j - 1] + sourceTerm) / (gammaUp + gammaDown);
                maxDelta = std::max(maxDelta, std::fabs(newK - k_[j]));
                k_[j] = std::max(newK, kFloor);
            }
        }

        // epsilon: interior cells solved; wall cells fixed at the
        // equilibrium wall-function value. The quadratic destruction term
        // (-Ceps2*epsilon^2/k) is linearized *implicitly* (Patankar's
        // positivity rule: sink terms go in the denominator, using the
        // previous sweep's epsilon as the frozen factor, not the numerator)
        // -- treating it fully explicitly caused runaway collapse during
        // development (a documented, common instability for this kind of
        // stiff quadratic sink, not specific to this implementation). The
        // production-like term stays an explicit numerator contribution,
        // per the same rule.
        for (int inner = 0; inner < 20; ++inner) {
            for (std::size_t j = 1; j + 1 < ny_; ++j) {
                const double gammaHere = nu_ + nut_[j] / kSigmaEps;
                const double gammaUp = 0.5 * (gammaHere + (nu_ + nut_[j + 1] / kSigmaEps));
                const double gammaDown = 0.5 * (gammaHere + (nu_ + nut_[j - 1] / kSigmaEps));
                const double kSafe = std::max(k_[j], kFloor);
                const double numerator = gammaUp * epsilon_[j + 1] + gammaDown * epsilon_[j - 1] +
                                          h_ * h_ * kCeps1 * (epsilon_[j] / kSafe) * production_[j];
                const double denominator =
                    (gammaUp + gammaDown) + h_ * h_ * kCeps2 * (epsilon_[j] / kSafe);
                const double newEps = numerator / denominator;
                maxDelta = std::max(maxDelta, std::fabs(newEps - epsilon_[j]));
                epsilon_[j] = std::max(newEps, kFloor);
            }
        }

        if (maxDelta < tolerance) {
            ++outer;
            break;
        }
    }
    return outer;
}

} // namespace aether::solver
