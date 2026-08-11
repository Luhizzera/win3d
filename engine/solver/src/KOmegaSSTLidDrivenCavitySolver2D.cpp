#include "aether/solver/KOmegaSSTLidDrivenCavitySolver2D.hpp"

#include "aether/solver/MixingLengthLidDrivenCavitySolver2D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

namespace {
constexpr double kKarman = 0.41;
constexpr double kBetaStar = 0.09;
constexpr double kA1 = 0.31;

constexpr double kSigmaK1 = 0.85;
constexpr double kSigmaW1 = 0.5;
constexpr double kBeta1 = 0.075;
constexpr double kSigmaK2 = 1.0;
constexpr double kSigmaW2 = 0.856;
constexpr double kBeta2 = 0.0828;

const double kGamma1 = kBeta1 / kBetaStar - kSigmaW1 * kKarman * kKarman / std::sqrt(kBetaStar);
const double kGamma2 = kBeta2 / kBetaStar - kSigmaW2 * kKarman * kKarman / std::sqrt(kBetaStar);

constexpr double kFloor = 1e-10;
} // namespace

KOmegaSSTLidDrivenCavitySolver2D::KOmegaSSTLidDrivenCavitySolver2D(std::size_t nx, std::size_t ny,
                                                                    double lengthX, double lengthY,
                                                                    double viscosity, double lidVelocity)
    : nx_(nx), ny_(ny), lengthX_(lengthX), lengthY_(lengthY), dx_(lengthX / static_cast<double>(nx)),
      dy_(lengthY / static_cast<double>(ny)), viscosity_(viscosity), lidVelocity_(lidVelocity),
      u_(nx * ny, 0.0), v_(nx * ny, 0.0), p_(nx * ny, 0.0), k_(nx * ny, 0.0), omega_(nx * ny, 0.0),
      nut_(nx * ny, 0.0), production_(nx * ny, 0.0), crossDiffusion_(nx * ny, 0.0),
      sigmaK_(nx * ny, kSigmaK2), sigmaOmega_(nx * ny, kSigmaW2), beta_(nx * ny, kBeta2),
      gamma_(nx * ny, kGamma2) {
    MixingLengthLidDrivenCavitySolver2D primer(nx, ny, lengthX, lengthY, viscosity, lidVelocity);
    const double primerDt = 0.3 * primer.stableTimeStep();
    for (int s = 0; s < 400; ++s) {
        primer.step(primerDt);
    }
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            u_[index(i, j)] = primer.u(i, j);
            v_[index(i, j)] = primer.v(i, j);
        }
    }

    const double k0 = 0.01 * lidVelocity * lidVelocity;
    const double lengthScale = 0.1 * std::min(lengthX_, lengthY_);
    // omega0 = eps0/(beta**k0) with eps0 = beta*^0.75*k0^1.5/L algebraically
    // simplifies to beta*^-0.25*sqrt(k0)/L -- written this way specifically
    // to avoid dividing by k0, which is exactly 0 when lidVelocity is 0
    // (the stationary-lid rest-state case): the naive eps0/(beta**k0) form
    // is a real 0/0 -> NaN this class hit during development, caught by
    // its own rest-state test failing with a NaN velocity field on the
    // very first step (nut_ picks up the NaN from omega_ before any
    // momentum diffusion term gets a chance to multiply it by an
    // identically-zero velocity difference).
    const double omega0 = std::sqrt(k0) / (lengthScale * std::pow(kBetaStar, 0.25));
    std::fill(k_.begin(), k_.end(), k0);
    std::fill(omega_.begin(), omega_.end(), omega0);
}

double KOmegaSSTLidDrivenCavitySolver2D::u(std::size_t i, std::size_t j) const { return u_[index(i, j)]; }
double KOmegaSSTLidDrivenCavitySolver2D::v(std::size_t i, std::size_t j) const { return v_[index(i, j)]; }
double KOmegaSSTLidDrivenCavitySolver2D::pressure(std::size_t i, std::size_t j) const {
    return p_[index(i, j)];
}
double KOmegaSSTLidDrivenCavitySolver2D::k(std::size_t i, std::size_t j) const { return k_[index(i, j)]; }
double KOmegaSSTLidDrivenCavitySolver2D::omega(std::size_t i, std::size_t j) const {
    return omega_[index(i, j)];
}
double KOmegaSSTLidDrivenCavitySolver2D::eddyViscosity(std::size_t i, std::size_t j) const {
    return nut_[index(i, j)];
}

double KOmegaSSTLidDrivenCavitySolver2D::lidVelocityAt(std::size_t i) const {
    constexpr double kPi = 3.14159265358979323846;
    const double x = (static_cast<double>(i) + 0.5) * dx_;
    const double s = std::sin(kPi * x / lengthX_);
    return lidVelocity_ * s * s;
}

double KOmegaSSTLidDrivenCavitySolver2D::wallDistanceAt(std::size_t i, std::size_t j) const {
    const double x = (static_cast<double>(i) + 0.5) * dx_;
    const double y = (static_cast<double>(j) + 0.5) * dy_;
    return std::min({x, lengthX_ - x, y, lengthY_ - y});
}

double KOmegaSSTLidDrivenCavitySolver2D::stableTimeStep() const {
    double maxNut = 0.0;
    for (double n : nut_) {
        maxNut = std::max(maxNut, n);
    }
    const double effectiveViscosity = viscosity_ + maxNut;
    const double diffusiveLimit =
        1.0 / (2.0 * effectiveViscosity * (1.0 / (dx_ * dx_) + 1.0 / (dy_ * dy_)));
    const double convectiveLimit = std::min(dx_, dy_) / std::max(lidVelocity_, 1e-12);
    return std::min(diffusiveLimit, convectiveLimit);
}

double KOmegaSSTLidDrivenCavitySolver2D::dirichletAt(const std::vector<double>& field, std::size_t i,
                                                      std::size_t j, int di, int dj,
                                                      double wallValue) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0) {
        return -field[index(0, j)];
    }
    if (ni >= static_cast<long long>(nx_)) {
        return -field[index(nx_ - 1, j)];
    }
    if (nj < 0) {
        return -field[index(i, 0)];
    }
    if (nj >= static_cast<long long>(ny_)) {
        return 2.0 * wallValue - field[index(i, ny_ - 1)];
    }
    return field[index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
}

double KOmegaSSTLidDrivenCavitySolver2D::neumannAt(const std::vector<double>& field, std::size_t i,
                                                    std::size_t j, int di, int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0) {
        return field[index(0, j)];
    }
    if (ni >= static_cast<long long>(nx_)) {
        return field[index(nx_ - 1, j)];
    }
    if (nj < 0) {
        return field[index(i, 0)];
    }
    if (nj >= static_cast<long long>(ny_)) {
        return field[index(i, ny_ - 1)];
    }
    return field[index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
}

double KOmegaSSTLidDrivenCavitySolver2D::omegaGhostAt(std::size_t i, std::size_t j, int di, int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    const bool crossesWall =
        ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 || nj >= static_cast<long long>(ny_);
    if (!crossesWall) {
        return omega_[index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
    }
    const double halfSpacing = (di != 0 ? dx_ : dy_) / 2.0;
    const double omegaWall = 2.0 * viscosity_ / (kBetaStar * halfSpacing * halfSpacing);
    return 2.0 * omegaWall - omega_[index(i, j)];
}

double KOmegaSSTLidDrivenCavitySolver2D::nutAt(std::size_t i, std::size_t j, int di, int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 || nj >= static_cast<long long>(ny_)) {
        return 0.0;
    }
    return nut_[index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
}

double KOmegaSSTLidDrivenCavitySolver2D::sigmaKNutAt(std::size_t i, std::size_t j, int di, int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 || nj >= static_cast<long long>(ny_)) {
        return 0.0;
    }
    const std::size_t idx = index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj));
    return sigmaK_[idx] * nut_[idx];
}

double KOmegaSSTLidDrivenCavitySolver2D::sigmaOmegaNutAt(std::size_t i, std::size_t j, int di,
                                                          int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 || nj >= static_cast<long long>(ny_)) {
        return 0.0;
    }
    const std::size_t idx = index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj));
    return sigmaOmega_[idx] * nut_[idx];
}

void KOmegaSSTLidDrivenCavitySolver2D::updateBlendingAndCoefficients() {
    constexpr double kRelaxation = 0.3;

    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            const double wallDist = std::max(wallDistanceAt(i, j), kFloor);
            const double kSafe = std::max(k_[idx], kFloor);
            const double omegaSafe = std::max(omega_[idx], kFloor);

            const double kE = dirichletAt(k_, i, j, 1, 0, 0.0);
            const double kW = dirichletAt(k_, i, j, -1, 0, 0.0);
            const double kN = dirichletAt(k_, i, j, 0, 1, 0.0);
            const double kS = dirichletAt(k_, i, j, 0, -1, 0.0);
            const double dkdx = (kE - kW) / (2.0 * dx_);
            const double dkdy = (kN - kS) / (2.0 * dy_);

            const double omegaE = omegaGhostAt(i, j, 1, 0);
            const double omegaW = omegaGhostAt(i, j, -1, 0);
            const double omegaN = omegaGhostAt(i, j, 0, 1);
            const double omegaS = omegaGhostAt(i, j, 0, -1);
            const double domegadx = (omegaE - omegaW) / (2.0 * dx_);
            const double domegady = (omegaN - omegaS) / (2.0 * dy_);

            const double gradDot = dkdx * domegadx + dkdy * domegady;
            const double cdKOmega = std::max(2.0 * kSigmaW2 / omegaSafe * gradDot, 1e-10);
            const double arg1a = std::sqrt(kSafe) / (kBetaStar * omegaSafe * wallDist);
            const double arg1b = 500.0 * viscosity_ / (wallDist * wallDist * omegaSafe);
            const double arg1c = 4.0 * kSigmaW2 * kSafe / (cdKOmega * wallDist * wallDist);
            const double arg1 = std::min(std::max(arg1a, arg1b), arg1c);
            const double f1 = std::tanh(arg1 * arg1 * arg1 * arg1);

            const double arg2a = 2.0 * std::sqrt(kSafe) / (kBetaStar * omegaSafe * wallDist);
            const double arg2b = 500.0 * viscosity_ / (wallDist * wallDist * omegaSafe);
            const double arg2 = std::max(arg2a, arg2b);
            const double f2 = std::tanh(arg2 * arg2);

            sigmaK_[idx] = f1 * kSigmaK1 + (1.0 - f1) * kSigmaK2;
            sigmaOmega_[idx] = f1 * kSigmaW1 + (1.0 - f1) * kSigmaW2;
            beta_[idx] = f1 * kBeta1 + (1.0 - f1) * kBeta2;
            gamma_[idx] = f1 * kGamma1 + (1.0 - f1) * kGamma2;

            const double uE = dirichletAt(u_, i, j, 1, 0, lidVelocityAt(i));
            const double uW = dirichletAt(u_, i, j, -1, 0, lidVelocityAt(i));
            const double uN = dirichletAt(u_, i, j, 0, 1, lidVelocityAt(i));
            const double uS = dirichletAt(u_, i, j, 0, -1, lidVelocityAt(i));
            const double vE = dirichletAt(v_, i, j, 1, 0, 0.0);
            const double vW = dirichletAt(v_, i, j, -1, 0, 0.0);
            const double vN = dirichletAt(v_, i, j, 0, 1, 0.0);
            const double vS = dirichletAt(v_, i, j, 0, -1, 0.0);

            const double dudx = (uE - uW) / (2.0 * dx_);
            const double dudy = (uN - uS) / (2.0 * dy_);
            const double dvdx = (vE - vW) / (2.0 * dx_);
            const double dvdy = (vN - vS) / (2.0 * dy_);

            const double shear = dudy + dvdx;
            const double strainSquared = 2.0 * dudx * dudx + 2.0 * dvdy * dvdy + shear * shear;
            // 2D vorticity magnitude (z-component of curl), used only for
            // the Bradshaw shear-stress limiter (Menter's own formula uses
            // vorticity there, not strain) -- distinct from the strain
            // rate used for production, though they coincide in the 1D
            // channel case where this class's 1D counterpart couldn't tell
            // them apart.
            const double vorticity = std::fabs(dvdx - dudy);

            const double nutComputed = kA1 * kSafe / std::max(kA1 * omegaSafe, vorticity * f2);
            nut_[idx] = kRelaxation * nutComputed + (1.0 - kRelaxation) * nut_[idx];

            const double productionRaw = nut_[idx] * strainSquared;
            production_[idx] = std::min(productionRaw, 10.0 * kBetaStar * kSafe * omegaSafe);

            crossDiffusion_[idx] = 2.0 * (1.0 - f1) * kSigmaW2 / omegaSafe * gradDot;
        }
    }
}

double KOmegaSSTLidDrivenCavitySolver2D::maxDivergence() const {
    double maxDiv = 0.0;
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const double uE = dirichletAt(u_, i, j, 1, 0, lidVelocityAt(i));
            const double uW = dirichletAt(u_, i, j, -1, 0, lidVelocityAt(i));
            const double vN = dirichletAt(v_, i, j, 0, 1, 0.0);
            const double vS = dirichletAt(v_, i, j, 0, -1, 0.0);
            const double div = (uE - uW) / (2.0 * dx_) + (vN - vS) / (2.0 * dy_);
            maxDiv = std::max(maxDiv, std::fabs(div));
        }
    }
    return maxDiv;
}

std::vector<double> KOmegaSSTLidDrivenCavitySolver2D::applyLaplacian(const std::vector<double>& x) const {
    const double ax = 1.0 / (dx_ * dx_);
    const double ay = 1.0 / (dy_ * dy_);
    const double weightTotal = 2.0 * ax + 2.0 * ay;

    std::vector<double> result(x.size());
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            if (idx == 0) {
                result[idx] = x[idx];
                continue;
            }
            const double xE = neumannAt(x, i, j, 1, 0);
            const double xW = neumannAt(x, i, j, -1, 0);
            const double xN = neumannAt(x, i, j, 0, 1);
            const double xS = neumannAt(x, i, j, 0, -1);
            const double weightedSum = ax * (xE + xW) + ay * (xN + xS);
            result[idx] = weightTotal * x[idx] - weightedSum;
        }
    }
    return result;
}

double KOmegaSSTLidDrivenCavitySolver2D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void KOmegaSSTLidDrivenCavitySolver2D::projectToDivergenceFree(std::vector<double>& uStar,
                                                                std::vector<double>& vStar, double dt) {
    const std::size_t n = uStar.size();

    std::vector<double> rhs(n);
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            if (idx == 0) {
                rhs[idx] = 0.0;
                continue;
            }
            const double uE = dirichletAt(uStar, i, j, 1, 0, lidVelocityAt(i));
            const double uW = dirichletAt(uStar, i, j, -1, 0, lidVelocityAt(i));
            const double vN = dirichletAt(vStar, i, j, 0, 1, 0.0);
            const double vS = dirichletAt(vStar, i, j, 0, -1, 0.0);
            const double divergence = (uE - uW) / (2.0 * dx_) + (vN - vS) / (2.0 * dy_);
            rhs[idx] = -divergence / dt;
        }
    }

    std::vector<double> residual(n);
    {
        const std::vector<double> operatorAppliedToP = applyLaplacian(p_);
        for (std::size_t idx = 0; idx < n; ++idx) {
            residual[idx] = idx == 0 ? 0.0 : rhs[idx] - operatorAppliedToP[idx];
        }
    }
    std::vector<double> direction = residual;
    double residualDotResidual = dot(residual, residual);

    for (std::size_t iteration = 0; iteration < n; ++iteration) {
        if (std::sqrt(residualDotResidual) < 1e-10) {
            break;
        }
        const std::vector<double> operatorAppliedToDirection = applyLaplacian(direction);
        const double directionDotOperatorDirection = dot(direction, operatorAppliedToDirection);
        if (directionDotOperatorDirection == 0.0) {
            break;
        }
        const double alpha = residualDotResidual / directionDotOperatorDirection;
        for (std::size_t idx = 0; idx < n; ++idx) {
            if (idx != 0) {
                p_[idx] += alpha * direction[idx];
            }
        }
        std::vector<double> newResidual(n);
        for (std::size_t idx = 0; idx < n; ++idx) {
            newResidual[idx] = idx == 0 ? 0.0 : residual[idx] - alpha * operatorAppliedToDirection[idx];
        }
        const double newResidualDotResidual = dot(newResidual, newResidual);
        const double beta = newResidualDotResidual / residualDotResidual;
        for (std::size_t idx = 0; idx < n; ++idx) {
            direction[idx] = idx == 0 ? 0.0 : newResidual[idx] + beta * direction[idx];
        }
        residual = std::move(newResidual);
        residualDotResidual = newResidualDotResidual;
    }

    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            const double pE = neumannAt(p_, i, j, 1, 0);
            const double pW = neumannAt(p_, i, j, -1, 0);
            const double pN = neumannAt(p_, i, j, 0, 1);
            const double pS = neumannAt(p_, i, j, 0, -1);
            const double gdpdx = (pE - pW) / (2.0 * dx_);
            const double gdpdy = (pN - pS) / (2.0 * dy_);
            u_[idx] = uStar[idx] - dt * gdpdx;
            v_[idx] = vStar[idx] - dt * gdpdy;
        }
    }
}

void KOmegaSSTLidDrivenCavitySolver2D::step(double dt) {
    updateBlendingAndCoefficients();

    std::vector<double> uStar(u_.size());
    std::vector<double> vStar(v_.size());
    std::vector<double> kStar(k_.size());
    std::vector<double> omegaStar(omega_.size());

    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);

            const double uE = dirichletAt(u_, i, j, 1, 0, lidVelocityAt(i));
            const double uW = dirichletAt(u_, i, j, -1, 0, lidVelocityAt(i));
            const double uN = dirichletAt(u_, i, j, 0, 1, lidVelocityAt(i));
            const double uS = dirichletAt(u_, i, j, 0, -1, lidVelocityAt(i));
            const double vE = dirichletAt(v_, i, j, 1, 0, 0.0);
            const double vW = dirichletAt(v_, i, j, -1, 0, 0.0);
            const double vN = dirichletAt(v_, i, j, 0, 1, 0.0);
            const double vS = dirichletAt(v_, i, j, 0, -1, 0.0);

            const double dudx = (uE - uW) / (2.0 * dx_);
            const double dudy = (uN - uS) / (2.0 * dy_);
            const double dvdx = (vE - vW) / (2.0 * dx_);
            const double dvdy = (vN - vS) / (2.0 * dy_);

            const double uConvection = u_[idx] * dudx + v_[idx] * dudy;
            const double vConvection = u_[idx] * dvdx + v_[idx] * dvdy;

            const double nutHere = nut_[idx];
            const double gammaE = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 1, 0));
            const double gammaW = viscosity_ + 0.5 * (nutHere + nutAt(i, j, -1, 0));
            const double gammaN = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 0, 1));
            const double gammaS = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 0, -1));

            const double diffusionU = (gammaE * (uE - u_[idx]) - gammaW * (u_[idx] - uW)) / (dx_ * dx_) +
                                       (gammaN * (uN - u_[idx]) - gammaS * (u_[idx] - uS)) / (dy_ * dy_);
            const double diffusionV = (gammaE * (vE - v_[idx]) - gammaW * (v_[idx] - vW)) / (dx_ * dx_) +
                                       (gammaN * (vN - v_[idx]) - gammaS * (v_[idx] - vS)) / (dy_ * dy_);

            uStar[idx] = u_[idx] + dt * (-uConvection + diffusionU);
            vStar[idx] = v_[idx] + dt * (-vConvection + diffusionV);

            // --- k transport: destruction beta**k*omega, Patankar-
            // linearized using the previous sweep's omega as the frozen
            // factor (same pattern as KEpsilonLidDrivenCavitySolver2D's
            // epsilon destruction). ---
            {
                const double kE = dirichletAt(k_, i, j, 1, 0, 0.0);
                const double kW = dirichletAt(k_, i, j, -1, 0, 0.0);
                const double kN = dirichletAt(k_, i, j, 0, 1, 0.0);
                const double kS = dirichletAt(k_, i, j, 0, -1, 0.0);
                const double dkdx = (kE - kW) / (2.0 * dx_);
                const double dkdy = (kN - kS) / (2.0 * dy_);
                const double advectionK = u_[idx] * dkdx + v_[idx] * dkdy;

                const double gEk = viscosity_ + 0.5 * (sigmaK_[idx] * nutHere + sigmaKNutAt(i, j, 1, 0));
                const double gWk = viscosity_ + 0.5 * (sigmaK_[idx] * nutHere + sigmaKNutAt(i, j, -1, 0));
                const double gNk = viscosity_ + 0.5 * (sigmaK_[idx] * nutHere + sigmaKNutAt(i, j, 0, 1));
                const double gSk = viscosity_ + 0.5 * (sigmaK_[idx] * nutHere + sigmaKNutAt(i, j, 0, -1));
                const double diffusionK = (gEk * (kE - k_[idx]) - gWk * (k_[idx] - kW)) / (dx_ * dx_) +
                                           (gNk * (kN - k_[idx]) - gSk * (k_[idx] - kS)) / (dy_ * dy_);

                const double explicitK =
                    k_[idx] + dt * (-advectionK + production_[idx] + diffusionK);
                const double newK = explicitK / (1.0 + dt * kBetaStar * omega_[idx]);
                kStar[idx] = std::max(newK, kFloor);
            }

            // --- omega transport: destruction beta*omega^2, Patankar-
            // linearized the same way. Cross-diffusion is an explicit
            // source term frozen for this step (computed in
            // updateBlendingAndCoefficients()). ---
            {
                const double omegaE = omegaGhostAt(i, j, 1, 0);
                const double omegaW = omegaGhostAt(i, j, -1, 0);
                const double omegaN = omegaGhostAt(i, j, 0, 1);
                const double omegaS = omegaGhostAt(i, j, 0, -1);
                const double domegadx = (omegaE - omegaW) / (2.0 * dx_);
                const double domegady = (omegaN - omegaS) / (2.0 * dy_);
                const double advectionOmega = u_[idx] * domegadx + v_[idx] * domegady;

                const double gEw =
                    viscosity_ + 0.5 * (sigmaOmega_[idx] * nutHere + sigmaOmegaNutAt(i, j, 1, 0));
                const double gWw =
                    viscosity_ + 0.5 * (sigmaOmega_[idx] * nutHere + sigmaOmegaNutAt(i, j, -1, 0));
                const double gNw =
                    viscosity_ + 0.5 * (sigmaOmega_[idx] * nutHere + sigmaOmegaNutAt(i, j, 0, 1));
                const double gSw =
                    viscosity_ + 0.5 * (sigmaOmega_[idx] * nutHere + sigmaOmegaNutAt(i, j, 0, -1));
                const double diffusionOmega =
                    (gEw * (omegaE - omega_[idx]) - gWw * (omega_[idx] - omegaW)) / (dx_ * dx_) +
                    (gNw * (omegaN - omega_[idx]) - gSw * (omega_[idx] - omegaS)) / (dy_ * dy_);

                const double nutSafe = std::max(nutHere, kFloor);
                const double productionOmega = gamma_[idx] / nutSafe * production_[idx];

                const double explicitOmega = omega_[idx] + dt * (-advectionOmega + productionOmega +
                                                                    crossDiffusion_[idx] + diffusionOmega);
                const double newOmega = explicitOmega / (1.0 + dt * beta_[idx] * omega_[idx]);
                omegaStar[idx] = std::max(newOmega, kFloor);
            }
        }
    }

    k_ = std::move(kStar);
    omega_ = std::move(omegaStar);

    projectToDivergenceFree(uStar, vStar, dt);
    time_ += dt;
}

} // namespace aether::solver
