#include "aether/solver/KOmegaSSTChannelFlowSolver1D.hpp"

#include "aether/solver/MixingLengthChannelFlowSolver1D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

namespace {
constexpr double kKarman = 0.41;
constexpr double kBetaStar = 0.09; // same value as k-epsilon's C_mu
constexpr double kA1 = 0.31;       // Bradshaw shear-stress limiter constant

// Set 1 (inner, Wilcox k-omega) and set 2 (outer, k-epsilon-derived).
constexpr double kSigmaK1 = 0.85;
constexpr double kSigmaW1 = 0.5;
constexpr double kBeta1 = 0.075;
constexpr double kSigmaK2 = 1.0;
constexpr double kSigmaW2 = 0.856;
constexpr double kBeta2 = 0.0828;

// gamma_i = beta_i/beta* - sigma_wi*kappa^2/sqrt(beta*), standard SST relation.
const double kGamma1 = kBeta1 / kBetaStar - kSigmaW1 * kKarman * kKarman / std::sqrt(kBetaStar);
const double kGamma2 = kBeta2 / kBetaStar - kSigmaW2 * kKarman * kKarman / std::sqrt(kBetaStar);

constexpr double kLogLawB = 5.0; // standard log-law additive constant
constexpr double kFloor = 1e-12; // guards against division by ~0 during transient Picard iterates
} // namespace

KOmegaSSTChannelFlowSolver1D::KOmegaSSTChannelFlowSolver1D(std::size_t ny, double height,
                                                            double kinematicViscosity, double source)
    : ny_(ny), height_(height), h_(height / static_cast<double>(ny)), nu_(kinematicViscosity),
      source_(source), u_(ny, 0.0), k_(ny, 0.0), omega_(ny, 0.0), nut_(ny, 0.0), production_(ny, 0.0),
      crossDiffusion_(ny, 0.0), sigmaK_(ny, kSigmaK2), sigmaOmega_(ny, kSigmaW2), beta_(ny, kBeta2),
      gamma_(ny, kGamma2) {
    const double uTau = frictionVelocity();
    const double kWall = uTau * uTau / std::sqrt(kBetaStar);
    for (std::size_t j = 0; j < ny_; ++j) {
        k_[j] = kWall;
        const double y = std::max(wallDistance(j), h_ * 0.5);
        omega_[j] = uTau / (kKarman * y * std::sqrt(kBetaStar));
    }

    // Warm-start the velocity field -- same chicken-and-egg deadlock as
    // KEpsilonChannelFlowSolver1D: u=0 gives zero production, so k/omega
    // collapse before u ever develops.
    MixingLengthChannelFlowSolver1D primer(ny, height, kinematicViscosity, source);
    primer.solve();
    for (std::size_t j = 0; j < ny_; ++j) {
        u_[j] = primer.u(j);
    }
}

double KOmegaSSTChannelFlowSolver1D::cellCenterY(std::size_t j) const {
    return (static_cast<double>(j) + 0.5) * h_;
}

double KOmegaSSTChannelFlowSolver1D::wallDistance(std::size_t j) const {
    const double y = cellCenterY(j);
    return std::min(y, height_ - y);
}

double KOmegaSSTChannelFlowSolver1D::frictionVelocity() const {
    return std::sqrt(source_ * height_ / 2.0);
}

double KOmegaSSTChannelFlowSolver1D::velocityGradientAt(std::size_t j) const {
    return (u_[j + 1] - u_[j - 1]) / (2.0 * h_);
}

void KOmegaSSTChannelFlowSolver1D::updateBlendingAndCoefficients() {
    // Under-relaxed for the same reason as KEpsilonChannelFlowSolver1D's
    // nu_t: it feeds back into both production and the transport
    // equations' own diffusivity, and without relaxation this loop
    // oscillates rather than settling.
    constexpr double kRelaxation = 0.3;

    for (std::size_t j = 1; j + 1 < ny_; ++j) {
        const double y = std::max(wallDistance(j), kFloor);
        const double kSafe = std::max(k_[j], kFloor);
        const double omegaSafe = std::max(omega_[j], kFloor);
        const double dkdy = (k_[j + 1] - k_[j - 1]) / (2.0 * h_);
        const double domegady = (omega_[j + 1] - omega_[j - 1]) / (2.0 * h_);

        const double cdKOmega = std::max(2.0 * kSigmaW2 / omegaSafe * dkdy * domegady, 1e-10);
        const double arg1a = std::sqrt(kSafe) / (kBetaStar * omegaSafe * y);
        const double arg1b = 500.0 * nu_ / (y * y * omegaSafe);
        const double arg1c = 4.0 * kSigmaW2 * kSafe / (cdKOmega * y * y);
        const double arg1 = std::min(std::max(arg1a, arg1b), arg1c);
        const double f1 = std::tanh(arg1 * arg1 * arg1 * arg1);

        const double arg2a = 2.0 * std::sqrt(kSafe) / (kBetaStar * omegaSafe * y);
        const double arg2b = 500.0 * nu_ / (y * y * omegaSafe);
        const double arg2 = std::max(arg2a, arg2b);
        const double f2 = std::tanh(arg2 * arg2);

        sigmaK_[j] = f1 * kSigmaK1 + (1.0 - f1) * kSigmaK2;
        sigmaOmega_[j] = f1 * kSigmaW1 + (1.0 - f1) * kSigmaW2;
        beta_[j] = f1 * kBeta1 + (1.0 - f1) * kBeta2;
        gamma_[j] = f1 * kGamma1 + (1.0 - f1) * kGamma2;

        const double gradient = velocityGradientAt(j);
        const double omegaMagnitude = std::fabs(gradient); // 1D shear: strain = vorticity magnitude
        const double nutComputed = kA1 * kSafe / std::max(kA1 * omegaSafe, omegaMagnitude * f2);
        nut_[j] = kRelaxation * nutComputed + (1.0 - kRelaxation) * nut_[j];

        const double productionRaw = nut_[j] * gradient * gradient;
        production_[j] = std::min(productionRaw, 10.0 * kBetaStar * kSafe * omegaSafe);

        crossDiffusion_[j] = 2.0 * (1.0 - f1) * kSigmaW2 / omegaSafe * dkdy * domegady;
    }
}

std::size_t KOmegaSSTChannelFlowSolver1D::solve(std::size_t maxOuterIterations, double tolerance) {
    const double uTau = frictionVelocity();
    const double kWall = uTau * uTau / std::sqrt(kBetaStar);
    const std::size_t last = ny_ - 1;
    const double omegaWallBottom = uTau / (kKarman * wallDistance(0) * std::sqrt(kBetaStar));
    const double omegaWallTop = uTau / (kKarman * wallDistance(last) * std::sqrt(kBetaStar));

    k_[0] = kWall;
    k_[last] = kWall;
    omega_[0] = omegaWallBottom;
    omega_[last] = omegaWallTop;

    // nu_t at the wall-function cells: updateBlendingAndCoefficients() only
    // covers interior cells (it needs a velocity-gradient stencil that
    // doesn't exist at the boundary), so these two would otherwise never be
    // set past their zero-initialized default -- but the wall-adjacent
    // momentum flux (gammaDown/gammaUp at j=1 and j=last-1) reads exactly
    // these entries, so leaving them at zero silently drops the turbulent
    // contribution to that one face. Fixed to the equilibrium mixing-length
    // value nu_t = kappa*u_tau*y_wall (the same value KEpsilonChannelFlowSolver1D
    // gets from Cmu*k_wall^2/eps_wall, and MixingLengthChannelFlowSolver1D
    // gets from l_m^2*|du/dy| in the log-law region -- all three reduce to
    // this identical formula under the equilibrium/log-law assumption
    // already baked into k_wall/omega_wall themselves) rather than
    // evaluating the SST blend formula itself, which needs a velocity
    // gradient and wall distance the boundary cell doesn't cleanly have.
    nut_[0] = kKarman * uTau * wallDistance(0);
    nut_[last] = kKarman * uTau * wallDistance(last);

    // Same null-space reasoning as KEpsilonChannelFlowSolver1D: fixing u at
    // the log-law value anchors the momentum equation, which a wall-shear
    // flux condition at both walls cannot do on its own.
    const double yPlusBottom = wallDistance(0) * uTau / nu_;
    const double yPlusTop = wallDistance(last) * uTau / nu_;
    u_[0] = uTau * ((1.0 / kKarman) * std::log(yPlusBottom) + kLogLawB);
    u_[last] = uTau * ((1.0 / kKarman) * std::log(yPlusTop) + kLogLawB);

    std::size_t outer = 0;
    for (; outer < maxOuterIterations; ++outer) {
        updateBlendingAndCoefficients();

        double maxDelta = 0.0;

        // Momentum: interior cells only; wall cells stay fixed at the
        // log-law value set above.
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

        // k: destruction (beta**k*omega) linearized implicitly using the
        // previous sweep's omega as the frozen factor -- a genuine product
        // of two transported unknowns, the same stiffness character as
        // k-epsilon's epsilon^2/k term.
        for (int inner = 0; inner < 20; ++inner) {
            for (std::size_t j = 1; j + 1 < ny_; ++j) {
                const double gammaHere = nu_ + sigmaK_[j] * nut_[j];
                const double gammaUp = 0.5 * (gammaHere + (nu_ + sigmaK_[j + 1] * nut_[j + 1]));
                const double gammaDown = 0.5 * (gammaHere + (nu_ + sigmaK_[j - 1] * nut_[j - 1]));
                const double numerator =
                    gammaUp * k_[j + 1] + gammaDown * k_[j - 1] + h_ * h_ * production_[j];
                const double denominator = (gammaUp + gammaDown) + h_ * h_ * kBetaStar * omega_[j];
                const double newK = numerator / denominator;
                maxDelta = std::max(maxDelta, std::fabs(newK - k_[j]));
                k_[j] = std::max(newK, kFloor);
            }
        }

        // omega: destruction (beta*omega^2) linearized implicitly the same
        // way as k-epsilon's epsilon destruction (using the previous
        // sweep's omega as one frozen factor). Cross-diffusion is an
        // explicit source term, frozen for the duration of this outer
        // iteration (computed in updateBlendingAndCoefficients()).
        for (int inner = 0; inner < 20; ++inner) {
            for (std::size_t j = 1; j + 1 < ny_; ++j) {
                const double gammaHere = nu_ + sigmaOmega_[j] * nut_[j];
                const double gammaUp = 0.5 * (gammaHere + (nu_ + sigmaOmega_[j + 1] * nut_[j + 1]));
                const double gammaDown = 0.5 * (gammaHere + (nu_ + sigmaOmega_[j - 1] * nut_[j - 1]));
                const double nutSafe = std::max(nut_[j], kFloor);
                const double numerator = gammaUp * omega_[j + 1] + gammaDown * omega_[j - 1] +
                                          h_ * h_ * (gamma_[j] / nutSafe * production_[j] +
                                                      crossDiffusion_[j]);
                const double denominator = (gammaUp + gammaDown) + h_ * h_ * beta_[j] * omega_[j];
                const double newOmega = numerator / denominator;
                maxDelta = std::max(maxDelta, std::fabs(newOmega - omega_[j]));
                omega_[j] = std::max(newOmega, kFloor);
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
