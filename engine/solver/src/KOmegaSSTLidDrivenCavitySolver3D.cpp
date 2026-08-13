#include "aether/solver/KOmegaSSTLidDrivenCavitySolver3D.hpp"

#include "aether/solver/MixingLengthLidDrivenCavitySolver3D.hpp"

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

KOmegaSSTLidDrivenCavitySolver3D::KOmegaSSTLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny,
                                                                     std::size_t nz, double lengthX,
                                                                     double lengthY, double lengthZ,
                                                                     double viscosity, double lidVelocity)
    : StaggeredCavityBase3D(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity),
      k_(nx * ny * nz, 0.0), omega_(nx * ny * nz, 0.0), nut_(nx * ny * nz, 0.0),
      production_(nx * ny * nz, 0.0), crossDiffusion_(nx * ny * nz, 0.0), sigmaK_(nx * ny * nz, kSigmaK2),
      sigmaOmega_(nx * ny * nz, kSigmaW2), beta_(nx * ny * nz, kBeta2), gamma_(nx * ny * nz, kGamma2) {
    // The base's momentum predictor and dt limit read nu_t through this.
    setEddyViscosityField(&nut_);

    MixingLengthLidDrivenCavitySolver3D primer(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity);
    const double primerDt = 0.3 * primer.stableTimeStep();
    for (int s = 0; s < 400; ++s) {
        primer.step(primerDt);
    }
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i <= nx_; ++i) {
                u_[indexU(i, j, k)] = primer.u(i, j, k);
            }
        }
    }
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j <= ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                v_[indexV(i, j, k)] = primer.v(i, j, k);
            }
        }
    }
    for (std::size_t k = 0; k <= nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                w_[indexW(i, j, k)] = primer.w(i, j, k);
            }
        }
    }

    const double k0 = 0.01 * lidVelocity * lidVelocity;
    const double lengthScale = 0.1 * std::min({lengthX_, lengthY_, lengthZ_});
    // omega0 = eps0/(beta**k0) with eps0 = beta*^0.75*k0^1.5/L simplifies
    // algebraically to beta*^-0.25*sqrt(k0)/L -- avoids the 0/0 -> NaN the
    // 2D class hit dividing by k0 (exactly 0 at lidVelocity=0) before this
    // simplification was applied. See KOmegaSSTLidDrivenCavitySolver2D's
    // own comment for the full story.
    const double omega0 = std::sqrt(k0) / (lengthScale * std::pow(kBetaStar, 0.25));
    std::fill(k_.begin(), k_.end(), k0);
    std::fill(omega_.begin(), omega_.end(), omega0);
}

double KOmegaSSTLidDrivenCavitySolver3D::k(std::size_t i, std::size_t j, std::size_t k) const {
    return k_[indexP(i, j, k)];
}

double KOmegaSSTLidDrivenCavitySolver3D::omega(std::size_t i, std::size_t j, std::size_t k) const {
    return omega_[indexP(i, j, k)];
}

double KOmegaSSTLidDrivenCavitySolver3D::eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const {
    return nut_[indexP(i, j, k)];
}

double KOmegaSSTLidDrivenCavitySolver3D::kAt(long long i, long long j, long long k) const {
    if (i < 0) {
        return -k_[indexP(0, static_cast<std::size_t>(j), static_cast<std::size_t>(k))];
    }
    if (i >= static_cast<long long>(nx_)) {
        return -k_[indexP(nx_ - 1, static_cast<std::size_t>(j), static_cast<std::size_t>(k))];
    }
    if (j < 0) {
        return -k_[indexP(static_cast<std::size_t>(i), 0, static_cast<std::size_t>(k))];
    }
    if (j >= static_cast<long long>(ny_)) {
        return -k_[indexP(static_cast<std::size_t>(i), ny_ - 1, static_cast<std::size_t>(k))];
    }
    if (k < 0) {
        return -k_[indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), 0)];
    }
    if (k >= static_cast<long long>(nz_)) {
        return -k_[indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), nz_ - 1)];
    }
    return k_[indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k))];
}

double KOmegaSSTLidDrivenCavitySolver3D::omegaGhostAt(std::size_t i, std::size_t j, std::size_t k, int di,
                                                        int dj, int dk) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    const long long nk = static_cast<long long>(k) + dk;
    const bool crossesWall = ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 ||
                              nj >= static_cast<long long>(ny_) || nk < 0 || nk >= static_cast<long long>(nz_);
    if (!crossesWall) {
        return omega_[indexP(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj),
                              static_cast<std::size_t>(nk))];
    }
    const double halfSpacing = (di != 0 ? dx_ : (dj != 0 ? dy_ : dz_)) / 2.0;
    const double omegaWall = 2.0 * viscosity_ / (kBetaStar * halfSpacing * halfSpacing);
    return 2.0 * omegaWall - omega_[indexP(i, j, k)];
}

double KOmegaSSTLidDrivenCavitySolver3D::sigmaKNutAt(long long i, long long j, long long k) const {
    if (i < 0 || i >= static_cast<long long>(nx_) || j < 0 || j >= static_cast<long long>(ny_) || k < 0 ||
        k >= static_cast<long long>(nz_)) {
        return 0.0;
    }
    const std::size_t idx =
        indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k));
    return sigmaK_[idx] * nut_[idx];
}

double KOmegaSSTLidDrivenCavitySolver3D::sigmaOmegaNutAt(long long i, long long j, long long k) const {
    if (i < 0 || i >= static_cast<long long>(nx_) || j < 0 || j >= static_cast<long long>(ny_) || k < 0 ||
        k >= static_cast<long long>(nz_)) {
        return 0.0;
    }
    const std::size_t idx =
        indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k));
    return sigmaOmega_[idx] * nut_[idx];
}

void KOmegaSSTLidDrivenCavitySolver3D::updateBlendingAndCoefficients() {
    constexpr double kRelaxation = 0.3;

    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = indexP(i, j, k);
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double wallDist = std::max(wallDistanceAt(i, j, k), kFloor);
                const double kSafe = std::max(k_[idx], kFloor);
                const double omegaSafe = std::max(omega_[idx], kFloor);

                const double kE = kAt(li + 1, lj, lk);
                const double kW = kAt(li - 1, lj, lk);
                const double kN = kAt(li, lj + 1, lk);
                const double kS = kAt(li, lj - 1, lk);
                const double kF = kAt(li, lj, lk + 1);
                const double kB = kAt(li, lj, lk - 1);
                const double dkdx = (kE - kW) / (2.0 * dx_);
                const double dkdy = (kN - kS) / (2.0 * dy_);
                const double dkdz = (kF - kB) / (2.0 * dz_);

                const double omegaE = omegaGhostAt(i, j, k, 1, 0, 0);
                const double omegaW = omegaGhostAt(i, j, k, -1, 0, 0);
                const double omegaN = omegaGhostAt(i, j, k, 0, 1, 0);
                const double omegaS = omegaGhostAt(i, j, k, 0, -1, 0);
                const double omegaF = omegaGhostAt(i, j, k, 0, 0, 1);
                const double omegaB = omegaGhostAt(i, j, k, 0, 0, -1);
                const double domegadx = (omegaE - omegaW) / (2.0 * dx_);
                const double domegady = (omegaN - omegaS) / (2.0 * dy_);
                const double domegadz = (omegaF - omegaB) / (2.0 * dz_);

                const double gradDot = dkdx * domegadx + dkdy * domegady + dkdz * domegadz;
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

                // Same 9-gradient block MixingLengthLidDrivenCavitySolver3D
                // uses for the strain-rate tensor, reused here for both
                // strain (production) and vorticity (Bradshaw limiter).
                const double dudx = (uAt(li + 1, lj, lk) - uAt(li, lj, lk)) / dx_;
                const double dvdy = (vAt(li, lj + 1, lk) - vAt(li, lj, lk)) / dy_;
                const double dwdz = (wAt(li, lj, lk + 1) - wAt(li, lj, lk)) / dz_;

                const double dudyFaceI = (uAt(li, lj + 1, lk) - uAt(li, lj - 1, lk)) / (2.0 * dy_);
                const double dudyFaceIp1 = (uAt(li + 1, lj + 1, lk) - uAt(li + 1, lj - 1, lk)) / (2.0 * dy_);
                const double dudy = 0.5 * (dudyFaceI + dudyFaceIp1);

                const double dvdxFaceJ = (vAt(li + 1, lj, lk) - vAt(li - 1, lj, lk)) / (2.0 * dx_);
                const double dvdxFaceJp1 = (vAt(li + 1, lj + 1, lk) - vAt(li - 1, lj + 1, lk)) / (2.0 * dx_);
                const double dvdx = 0.5 * (dvdxFaceJ + dvdxFaceJp1);

                const double dudzFaceI = (uAt(li, lj, lk + 1) - uAt(li, lj, lk - 1)) / (2.0 * dz_);
                const double dudzFaceIp1 = (uAt(li + 1, lj, lk + 1) - uAt(li + 1, lj, lk - 1)) / (2.0 * dz_);
                const double dudz = 0.5 * (dudzFaceI + dudzFaceIp1);

                const double dwdxFaceK = (wAt(li + 1, lj, lk) - wAt(li - 1, lj, lk)) / (2.0 * dx_);
                const double dwdxFaceKp1 = (wAt(li + 1, lj, lk + 1) - wAt(li - 1, lj, lk + 1)) / (2.0 * dx_);
                const double dwdx = 0.5 * (dwdxFaceK + dwdxFaceKp1);

                const double dvdzFaceJ = (vAt(li, lj, lk + 1) - vAt(li, lj, lk - 1)) / (2.0 * dz_);
                const double dvdzFaceJp1 = (vAt(li, lj + 1, lk + 1) - vAt(li, lj + 1, lk - 1)) / (2.0 * dz_);
                const double dvdz = 0.5 * (dvdzFaceJ + dvdzFaceJp1);

                const double dwdyFaceK = (wAt(li, lj + 1, lk) - wAt(li, lj - 1, lk)) / (2.0 * dy_);
                const double dwdyFaceKp1 = (wAt(li, lj + 1, lk + 1) - wAt(li, lj - 1, lk + 1)) / (2.0 * dy_);
                const double dwdy = 0.5 * (dwdyFaceK + dwdyFaceKp1);

                const double sxy = 0.5 * (dudy + dvdx);
                const double sxz = 0.5 * (dudz + dwdx);
                const double syz = 0.5 * (dvdz + dwdy);
                const double strainSquared = 2.0 * (dudx * dudx + dvdy * dvdy + dwdz * dwdz) +
                                              4.0 * (sxy * sxy + sxz * sxz + syz * syz);

                // Full 3D vorticity magnitude (curl(u)) for the Bradshaw
                // limiter -- the 2D class only ever saw its z-component
                // (dv/dx - du/dy), the sole nonzero one in a genuinely 2D
                // flow.
                const double vorticityX = dwdy - dvdz;
                const double vorticityY = dudz - dwdx;
                const double vorticityZ = dvdx - dudy;
                const double vorticity = std::sqrt(vorticityX * vorticityX + vorticityY * vorticityY +
                                                    vorticityZ * vorticityZ);

                const double nutComputed = kA1 * kSafe / std::max(kA1 * omegaSafe, vorticity * f2);
                nut_[idx] = kRelaxation * nutComputed + (1.0 - kRelaxation) * nut_[idx];

                const double productionRaw = nut_[idx] * strainSquared;
                production_[idx] = std::min(productionRaw, 10.0 * kBetaStar * kSafe * omegaSafe);

                crossDiffusion_[idx] = 2.0 * (1.0 - f1) * kSigmaW2 / omegaSafe * gradDot;
            }
        }
    }
}

void KOmegaSSTLidDrivenCavitySolver3D::step(double dt) {
    updateBlendingAndCoefficients();

    std::vector<double> uStar = u_;
    std::vector<double> vStar = v_;
    std::vector<double> wStar = w_;
    std::vector<double> kStar(k_.size());
    std::vector<double> omegaStar(omega_.size());

    computeMomentumPredictor(uStar, vStar, wStar, dt);

    // --- k and omega transport: colocated 6-neighbor stencil, both fields
    // living at cell centers like pressure and nu_t. ---
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = indexP(i, j, k);
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double uCenter = 0.5 * (uAt(li, lj, lk) + uAt(li + 1, lj, lk));
                const double vCenter = 0.5 * (vAt(li, lj, lk) + vAt(li, lj + 1, lk));
                const double wCenter = 0.5 * (wAt(li, lj, lk) + wAt(li, lj, lk + 1));

                const double nutHere = nut_[idx];
                const double sigmaKHere = sigmaK_[idx];
                const double sigmaOmegaHere = sigmaOmega_[idx];

                // --- k transport: destruction beta**k*omega, Patankar-
                // linearized using the current omega as the frozen factor.
                {
                    const double kE = kAt(li + 1, lj, lk);
                    const double kW = kAt(li - 1, lj, lk);
                    const double kN = kAt(li, lj + 1, lk);
                    const double kS = kAt(li, lj - 1, lk);
                    const double kF = kAt(li, lj, lk + 1);
                    const double kB = kAt(li, lj, lk - 1);
                    const double kHere = k_[idx];

                    const double dkdx = (kE - kW) / (2.0 * dx_);
                    const double dkdy = (kN - kS) / (2.0 * dy_);
                    const double dkdz = (kF - kB) / (2.0 * dz_);
                    const double advectionK = uCenter * dkdx + vCenter * dkdy + wCenter * dkdz;

                    const double gEk = viscosity_ + 0.5 * (sigmaKHere * nutHere + sigmaKNutAt(li + 1, lj, lk));
                    const double gWk = viscosity_ + 0.5 * (sigmaKHere * nutHere + sigmaKNutAt(li - 1, lj, lk));
                    const double gNk = viscosity_ + 0.5 * (sigmaKHere * nutHere + sigmaKNutAt(li, lj + 1, lk));
                    const double gSk = viscosity_ + 0.5 * (sigmaKHere * nutHere + sigmaKNutAt(li, lj - 1, lk));
                    const double gFk = viscosity_ + 0.5 * (sigmaKHere * nutHere + sigmaKNutAt(li, lj, lk + 1));
                    const double gBk = viscosity_ + 0.5 * (sigmaKHere * nutHere + sigmaKNutAt(li, lj, lk - 1));

                    const double diffusionK = (gEk * (kE - kHere) - gWk * (kHere - kW)) / (dx_ * dx_) +
                                               (gNk * (kN - kHere) - gSk * (kHere - kS)) / (dy_ * dy_) +
                                               (gFk * (kF - kHere) - gBk * (kHere - kB)) / (dz_ * dz_);

                    const double explicitK = kHere + dt * (-advectionK + production_[idx] + diffusionK);
                    const double newK = explicitK / (1.0 + dt * kBetaStar * omega_[idx]);
                    kStar[idx] = std::max(newK, kFloor);
                }

                // --- omega transport: destruction beta*omega^2, Patankar-
                // linearized the same way; cross-diffusion is an explicit
                // frozen source term from updateBlendingAndCoefficients().
                {
                    const double omegaE = omegaGhostAt(i, j, k, 1, 0, 0);
                    const double omegaW = omegaGhostAt(i, j, k, -1, 0, 0);
                    const double omegaN = omegaGhostAt(i, j, k, 0, 1, 0);
                    const double omegaS = omegaGhostAt(i, j, k, 0, -1, 0);
                    const double omegaF = omegaGhostAt(i, j, k, 0, 0, 1);
                    const double omegaB = omegaGhostAt(i, j, k, 0, 0, -1);
                    const double omegaHere = omega_[idx];

                    const double domegadx = (omegaE - omegaW) / (2.0 * dx_);
                    const double domegady = (omegaN - omegaS) / (2.0 * dy_);
                    const double domegadz = (omegaF - omegaB) / (2.0 * dz_);
                    const double advectionOmega = uCenter * domegadx + vCenter * domegady + wCenter * domegadz;

                    const double gEw =
                        viscosity_ + 0.5 * (sigmaOmegaHere * nutHere + sigmaOmegaNutAt(li + 1, lj, lk));
                    const double gWw =
                        viscosity_ + 0.5 * (sigmaOmegaHere * nutHere + sigmaOmegaNutAt(li - 1, lj, lk));
                    const double gNw =
                        viscosity_ + 0.5 * (sigmaOmegaHere * nutHere + sigmaOmegaNutAt(li, lj + 1, lk));
                    const double gSw =
                        viscosity_ + 0.5 * (sigmaOmegaHere * nutHere + sigmaOmegaNutAt(li, lj - 1, lk));
                    const double gFw =
                        viscosity_ + 0.5 * (sigmaOmegaHere * nutHere + sigmaOmegaNutAt(li, lj, lk + 1));
                    const double gBw =
                        viscosity_ + 0.5 * (sigmaOmegaHere * nutHere + sigmaOmegaNutAt(li, lj, lk - 1));

                    const double diffusionOmega =
                        (gEw * (omegaE - omegaHere) - gWw * (omegaHere - omegaW)) / (dx_ * dx_) +
                        (gNw * (omegaN - omegaHere) - gSw * (omegaHere - omegaS)) / (dy_ * dy_) +
                        (gFw * (omegaF - omegaHere) - gBw * (omegaHere - omegaB)) / (dz_ * dz_);

                    const double nutSafe = std::max(nutHere, kFloor);
                    const double productionOmega = gamma_[idx] / nutSafe * production_[idx];

                    const double explicitOmega = omegaHere + dt * (-advectionOmega + productionOmega +
                                                                      crossDiffusion_[idx] + diffusionOmega);
                    const double newOmega = explicitOmega / (1.0 + dt * beta_[idx] * omegaHere);
                    omegaStar[idx] = std::max(newOmega, kFloor);
                }
            }
        }
    }

    k_ = std::move(kStar);
    omega_ = std::move(omegaStar);

    projectToDivergenceFree(uStar, vStar, wStar, dt);
    time_ += dt;
}

} // namespace aether::solver
