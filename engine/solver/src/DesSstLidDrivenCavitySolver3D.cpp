#include "aether/solver/DesSstLidDrivenCavitySolver3D.hpp"

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

DesSstLidDrivenCavitySolver3D::DesSstLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz,
                                                               double lengthX, double lengthY, double lengthZ,
                                                               double viscosity, double lidVelocity, double cDes)
    : nx_(nx), ny_(ny), nz_(nz), lengthX_(lengthX), lengthY_(lengthY), lengthZ_(lengthZ),
      dx_(lengthX / static_cast<double>(nx)), dy_(lengthY / static_cast<double>(ny)),
      dz_(lengthZ / static_cast<double>(nz)), viscosity_(viscosity), lidVelocity_(lidVelocity), cDes_(cDes),
      u_((nx + 1) * ny * nz, 0.0), v_(nx * (ny + 1) * nz, 0.0), w_(nx * ny * (nz + 1), 0.0),
      p_(nx * ny * nz, 0.0), k_(nx * ny * nz, 0.0), omega_(nx * ny * nz, 0.0), nut_(nx * ny * nz, 0.0),
      production_(nx * ny * nz, 0.0), crossDiffusion_(nx * ny * nz, 0.0), sigmaK_(nx * ny * nz, kSigmaK2),
      sigmaOmega_(nx * ny * nz, kSigmaW2), beta_(nx * ny * nz, kBeta2), gamma_(nx * ny * nz, kGamma2),
      fDes_(nx * ny * nz, 1.0) {
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
    const double omega0 = std::sqrt(k0) / (lengthScale * std::pow(kBetaStar, 0.25));
    std::fill(k_.begin(), k_.end(), k0);
    std::fill(omega_.begin(), omega_.end(), omega0);
}

double DesSstLidDrivenCavitySolver3D::u(std::size_t i, std::size_t j, std::size_t k) const {
    return u_[indexU(i, j, k)];
}
double DesSstLidDrivenCavitySolver3D::v(std::size_t i, std::size_t j, std::size_t k) const {
    return v_[indexV(i, j, k)];
}
double DesSstLidDrivenCavitySolver3D::w(std::size_t i, std::size_t j, std::size_t k) const {
    return w_[indexW(i, j, k)];
}
double DesSstLidDrivenCavitySolver3D::pressure(std::size_t i, std::size_t j, std::size_t k) const {
    return p_[indexP(i, j, k)];
}
double DesSstLidDrivenCavitySolver3D::k(std::size_t i, std::size_t j, std::size_t k) const {
    return k_[indexP(i, j, k)];
}
double DesSstLidDrivenCavitySolver3D::omega(std::size_t i, std::size_t j, std::size_t k) const {
    return omega_[indexP(i, j, k)];
}
double DesSstLidDrivenCavitySolver3D::eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const {
    return nut_[indexP(i, j, k)];
}
double DesSstLidDrivenCavitySolver3D::desFactor(std::size_t i, std::size_t j, std::size_t k) const {
    return fDes_[indexP(i, j, k)];
}
double DesSstLidDrivenCavitySolver3D::filterWidth() const { return std::cbrt(dx_ * dy_ * dz_); }

double DesSstLidDrivenCavitySolver3D::lidVelocityAt(double x, double y) const {
    constexpr double kPi = 3.14159265358979323846;
    const double sx = std::sin(kPi * x / lengthX_);
    const double sy = std::sin(kPi * y / lengthY_);
    return lidVelocity_ * sx * sx * sy * sy;
}

double DesSstLidDrivenCavitySolver3D::wallDistanceAt(std::size_t i, std::size_t j, std::size_t k) const {
    const double x = (static_cast<double>(i) + 0.5) * dx_;
    const double y = (static_cast<double>(j) + 0.5) * dy_;
    const double z = (static_cast<double>(k) + 0.5) * dz_;
    return std::min({x, lengthX_ - x, y, lengthY_ - y, z, lengthZ_ - z});
}

double DesSstLidDrivenCavitySolver3D::stableTimeStep() const {
    double maxNut = 0.0;
    for (double n : nut_) {
        maxNut = std::max(maxNut, n);
    }
    const double effectiveViscosity = viscosity_ + maxNut;
    const double diffusiveLimit = 1.0 / (2.0 * effectiveViscosity *
                                          (1.0 / (dx_ * dx_) + 1.0 / (dy_ * dy_) + 1.0 / (dz_ * dz_)));
    const double convectiveLimit = std::min({dx_, dy_, dz_}) / std::max(lidVelocity_, 1e-12);
    return std::min(diffusiveLimit, convectiveLimit);
}

double DesSstLidDrivenCavitySolver3D::uAt(long long i, long long j, long long k) const {
    const auto ii = static_cast<std::size_t>(i);
    if (j < 0) {
        return -u_[indexU(ii, 0, static_cast<std::size_t>(k))];
    }
    if (j >= static_cast<long long>(ny_)) {
        return -u_[indexU(ii, ny_ - 1, static_cast<std::size_t>(k))];
    }
    if (k < 0) {
        return -u_[indexU(ii, static_cast<std::size_t>(j), 0)];
    }
    if (k >= static_cast<long long>(nz_)) {
        const double x = static_cast<double>(i) * dx_;
        const double y = (static_cast<double>(j) + 0.5) * dy_;
        return 2.0 * lidVelocityAt(x, y) - u_[indexU(ii, static_cast<std::size_t>(j), nz_ - 1)];
    }
    return u_[indexU(ii, static_cast<std::size_t>(j), static_cast<std::size_t>(k))];
}

double DesSstLidDrivenCavitySolver3D::vAt(long long i, long long j, long long k) const {
    const auto jj = static_cast<std::size_t>(j);
    if (i < 0) {
        return -v_[indexV(0, jj, static_cast<std::size_t>(k))];
    }
    if (i >= static_cast<long long>(nx_)) {
        return -v_[indexV(nx_ - 1, jj, static_cast<std::size_t>(k))];
    }
    if (k < 0) {
        return -v_[indexV(static_cast<std::size_t>(i), jj, 0)];
    }
    if (k >= static_cast<long long>(nz_)) {
        return -v_[indexV(static_cast<std::size_t>(i), jj, nz_ - 1)];
    }
    return v_[indexV(static_cast<std::size_t>(i), jj, static_cast<std::size_t>(k))];
}

double DesSstLidDrivenCavitySolver3D::wAt(long long i, long long j, long long k) const {
    const auto kk = static_cast<std::size_t>(k);
    if (i < 0) {
        return -w_[indexW(0, static_cast<std::size_t>(j), kk)];
    }
    if (i >= static_cast<long long>(nx_)) {
        return -w_[indexW(nx_ - 1, static_cast<std::size_t>(j), kk)];
    }
    if (j < 0) {
        return -w_[indexW(static_cast<std::size_t>(i), 0, kk)];
    }
    if (j >= static_cast<long long>(ny_)) {
        return -w_[indexW(static_cast<std::size_t>(i), ny_ - 1, kk)];
    }
    return w_[indexW(static_cast<std::size_t>(i), static_cast<std::size_t>(j), kk)];
}

double DesSstLidDrivenCavitySolver3D::pAt(long long i, long long j, long long k) const {
    long long ci = i;
    long long cj = j;
    long long ck = k;
    if (ci < 0) {
        ci = 0;
    } else if (ci >= static_cast<long long>(nx_)) {
        ci = static_cast<long long>(nx_) - 1;
    }
    if (cj < 0) {
        cj = 0;
    } else if (cj >= static_cast<long long>(ny_)) {
        cj = static_cast<long long>(ny_) - 1;
    }
    if (ck < 0) {
        ck = 0;
    } else if (ck >= static_cast<long long>(nz_)) {
        ck = static_cast<long long>(nz_) - 1;
    }
    return p_[indexP(static_cast<std::size_t>(ci), static_cast<std::size_t>(cj), static_cast<std::size_t>(ck))];
}

double DesSstLidDrivenCavitySolver3D::nutAt(long long i, long long j, long long k) const {
    if (i < 0 || i >= static_cast<long long>(nx_) || j < 0 || j >= static_cast<long long>(ny_) || k < 0 ||
        k >= static_cast<long long>(nz_)) {
        return 0.0;
    }
    return nut_[indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k))];
}

double DesSstLidDrivenCavitySolver3D::kAt(long long i, long long j, long long k) const {
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

double DesSstLidDrivenCavitySolver3D::omegaGhostAt(std::size_t i, std::size_t j, std::size_t k, int di, int dj,
                                                     int dk) const {
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

double DesSstLidDrivenCavitySolver3D::sigmaKNutAt(long long i, long long j, long long k) const {
    if (i < 0 || i >= static_cast<long long>(nx_) || j < 0 || j >= static_cast<long long>(ny_) || k < 0 ||
        k >= static_cast<long long>(nz_)) {
        return 0.0;
    }
    const std::size_t idx =
        indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k));
    return sigmaK_[idx] * nut_[idx];
}

double DesSstLidDrivenCavitySolver3D::sigmaOmegaNutAt(long long i, long long j, long long k) const {
    if (i < 0 || i >= static_cast<long long>(nx_) || j < 0 || j >= static_cast<long long>(ny_) || k < 0 ||
        k >= static_cast<long long>(nz_)) {
        return 0.0;
    }
    const std::size_t idx =
        indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k));
    return sigmaOmega_[idx] * nut_[idx];
}

void DesSstLidDrivenCavitySolver3D::updateBlendingAndCoefficients() {
    constexpr double kRelaxation = 0.3;
    const double delta = filterWidth();

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

                // DES length-scale substitution: L_RANS = sqrt(k)/(beta*_star*omega)
                // (the RANS turbulence length scale implicit in the ordinary
                // k-equation destruction term beta_star*k*omega = k^1.5/L_RANS).
                // F_DES = max(L_RANS/(C_DES*Delta), 1.0) multiplies that
                // destruction term in step(); F_DES=1 leaves plain SST
                // unchanged, F_DES>1 is the LES-mode switch. See the class
                // header for the full derivation.
                const double lRans = std::sqrt(kSafe) / (kBetaStar * omegaSafe);
                fDes_[idx] = std::max(lRans / (cDes_ * delta), 1.0);

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

double DesSstLidDrivenCavitySolver3D::maxDivergence() const {
    double maxDiv = 0.0;
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const double div = (u_[indexU(i + 1, j, k)] - u_[indexU(i, j, k)]) / dx_ +
                                   (v_[indexV(i, j + 1, k)] - v_[indexV(i, j, k)]) / dy_ +
                                   (w_[indexW(i, j, k + 1)] - w_[indexW(i, j, k)]) / dz_;
                maxDiv = std::max(maxDiv, std::fabs(div));
            }
        }
    }
    return maxDiv;
}

std::vector<double> DesSstLidDrivenCavitySolver3D::applyLaplacian(const std::vector<double>& x) const {
    const double ax = 1.0 / (dx_ * dx_);
    const double ay = 1.0 / (dy_ * dy_);
    const double az = 1.0 / (dz_ * dz_);
    const double weightTotal = 2.0 * ax + 2.0 * ay + 2.0 * az;

    std::vector<double> result(x.size());
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = indexP(i, j, k);
                if (idx == 0) {
                    result[idx] = x[idx];
                    continue;
                }
                const double left = i > 0 ? x[indexP(i - 1, j, k)] : x[indexP(0, j, k)];
                const double right = i + 1 < nx_ ? x[indexP(i + 1, j, k)] : x[indexP(nx_ - 1, j, k)];
                const double down = j > 0 ? x[indexP(i, j - 1, k)] : x[indexP(i, 0, k)];
                const double up = j + 1 < ny_ ? x[indexP(i, j + 1, k)] : x[indexP(i, ny_ - 1, k)];
                const double back = k > 0 ? x[indexP(i, j, k - 1)] : x[indexP(i, j, 0)];
                const double front = k + 1 < nz_ ? x[indexP(i, j, k + 1)] : x[indexP(i, j, nz_ - 1)];
                const double weightedSum = ax * (left + right) + ay * (down + up) + az * (back + front);
                result[idx] = weightTotal * x[idx] - weightedSum;
            }
        }
    }
    return result;
}

double DesSstLidDrivenCavitySolver3D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void DesSstLidDrivenCavitySolver3D::projectToDivergenceFree(std::vector<double>& uStar,
                                                              std::vector<double>& vStar,
                                                              std::vector<double>& wStar, double dt) {
    const std::size_t n = p_.size();

    std::vector<double> rhs(n);
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = indexP(i, j, k);
                if (idx == 0) {
                    rhs[idx] = 0.0;
                    continue;
                }
                const double divergence = (uStar[indexU(i + 1, j, k)] - uStar[indexU(i, j, k)]) / dx_ +
                                           (vStar[indexV(i, j + 1, k)] - vStar[indexV(i, j, k)]) / dy_ +
                                           (wStar[indexW(i, j, k + 1)] - wStar[indexW(i, j, k)]) / dz_;
                rhs[idx] = -divergence / dt;
            }
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

    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 1; i < nx_; ++i) {
                const double dpdx = (p_[indexP(i, j, k)] - p_[indexP(i - 1, j, k)]) / dx_;
                u_[indexU(i, j, k)] = uStar[indexU(i, j, k)] - dt * dpdx;
            }
        }
    }
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 1; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const double dpdy = (p_[indexP(i, j, k)] - p_[indexP(i, j - 1, k)]) / dy_;
                v_[indexV(i, j, k)] = vStar[indexV(i, j, k)] - dt * dpdy;
            }
        }
    }
    for (std::size_t k = 1; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const double dpdz = (p_[indexP(i, j, k)] - p_[indexP(i, j, k - 1)]) / dz_;
                w_[indexW(i, j, k)] = wStar[indexW(i, j, k)] - dt * dpdz;
            }
        }
    }
}

void DesSstLidDrivenCavitySolver3D::step(double dt) {
    updateBlendingAndCoefficients();

    std::vector<double> uStar = u_;
    std::vector<double> vStar = v_;
    std::vector<double> wStar = w_;
    std::vector<double> kStar(k_.size());
    std::vector<double> omegaStar(omega_.size());

    // --- u-momentum: interior x-faces only (i in [1, nx-1]) ---
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 1; i < nx_; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double uHere = uAt(li, lj, lk);
                const double uCenterAtI = 0.5 * (uHere + uAt(li + 1, lj, lk));
                const double uCenterAtIm1 = 0.5 * (uAt(li - 1, lj, lk) + uHere);
                const double duudx = (uCenterAtI * uCenterAtI - uCenterAtIm1 * uCenterAtIm1) / dx_;

                const double uEdgeJp = 0.5 * (uHere + uAt(li, lj + 1, lk));
                const double uEdgeJm = 0.5 * (uAt(li, lj - 1, lk) + uHere);
                const double vEdgeIJp = 0.5 * (vAt(li - 1, lj + 1, lk) + vAt(li, lj + 1, lk));
                const double vEdgeIJ = 0.5 * (vAt(li - 1, lj, lk) + vAt(li, lj, lk));
                const double duvdy = (uEdgeJp * vEdgeIJp - uEdgeJm * vEdgeIJ) / dy_;

                const double uEdgeKp = 0.5 * (uHere + uAt(li, lj, lk + 1));
                const double uEdgeKm = 0.5 * (uAt(li, lj, lk - 1) + uHere);
                const double wEdgeIKp = 0.5 * (wAt(li - 1, lj, lk + 1) + wAt(li, lj, lk + 1));
                const double wEdgeIK = 0.5 * (wAt(li - 1, lj, lk) + wAt(li, lj, lk));
                const double duwdz = (uEdgeKp * wEdgeIKp - uEdgeKm * wEdgeIK) / dz_;

                const double gammaE = viscosity_ + nutAt(li, lj, lk);
                const double gammaW = viscosity_ + nutAt(li - 1, lj, lk);
                const double gammaTransverse = viscosity_ + 0.5 * (nutAt(li - 1, lj, lk) + nutAt(li, lj, lk));

                const double diffusionU =
                    (gammaE * (uAt(li + 1, lj, lk) - uHere) - gammaW * (uHere - uAt(li - 1, lj, lk))) /
                        (dx_ * dx_) +
                    gammaTransverse * (uAt(li, lj + 1, lk) - 2.0 * uHere + uAt(li, lj - 1, lk)) / (dy_ * dy_) +
                    gammaTransverse * (uAt(li, lj, lk + 1) - 2.0 * uHere + uAt(li, lj, lk - 1)) / (dz_ * dz_);

                uStar[indexU(i, j, k)] = uHere + dt * (-(duudx + duvdy + duwdz) + diffusionU);
            }
        }
    }

    // --- v-momentum: interior y-faces only (j in [1, ny-1]) ---
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 1; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double vHere = vAt(li, lj, lk);
                const double vCenterAtJ = 0.5 * (vHere + vAt(li, lj + 1, lk));
                const double vCenterAtJm1 = 0.5 * (vAt(li, lj - 1, lk) + vHere);
                const double dvvdy = (vCenterAtJ * vCenterAtJ - vCenterAtJm1 * vCenterAtJm1) / dy_;

                const double vEdgeKp = 0.5 * (vHere + vAt(li, lj, lk + 1));
                const double vEdgeKm = 0.5 * (vAt(li, lj, lk - 1) + vHere);
                const double wEdgeJKp = 0.5 * (wAt(li, lj - 1, lk + 1) + wAt(li, lj, lk + 1));
                const double wEdgeJK = 0.5 * (wAt(li, lj - 1, lk) + wAt(li, lj, lk));
                const double dvwdz = (vEdgeKp * wEdgeJKp - vEdgeKm * wEdgeJK) / dz_;

                const double vEdgeIp = 0.5 * (vHere + vAt(li + 1, lj, lk));
                const double vEdgeIm = 0.5 * (vAt(li - 1, lj, lk) + vHere);
                const double uEdgeJIp = 0.5 * (uAt(li + 1, lj - 1, lk) + uAt(li + 1, lj, lk));
                const double uEdgeJI = 0.5 * (uAt(li, lj - 1, lk) + uAt(li, lj, lk));
                const double dvudx = (vEdgeIp * uEdgeJIp - vEdgeIm * uEdgeJI) / dx_;

                const double gammaN = viscosity_ + nutAt(li, lj, lk);
                const double gammaS = viscosity_ + nutAt(li, lj - 1, lk);
                const double gammaTransverse = viscosity_ + 0.5 * (nutAt(li, lj - 1, lk) + nutAt(li, lj, lk));

                const double diffusionV =
                    (gammaN * (vAt(li, lj + 1, lk) - vHere) - gammaS * (vHere - vAt(li, lj - 1, lk))) /
                        (dy_ * dy_) +
                    gammaTransverse * (vAt(li + 1, lj, lk) - 2.0 * vHere + vAt(li - 1, lj, lk)) / (dx_ * dx_) +
                    gammaTransverse * (vAt(li, lj, lk + 1) - 2.0 * vHere + vAt(li, lj, lk - 1)) / (dz_ * dz_);

                vStar[indexV(i, j, k)] = vHere + dt * (-(dvvdy + dvwdz + dvudx) + diffusionV);
            }
        }
    }

    // --- w-momentum: interior z-faces only (k in [1, nz-1]) ---
    for (std::size_t k = 1; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double wHere = wAt(li, lj, lk);
                const double wCenterAtK = 0.5 * (wHere + wAt(li, lj, lk + 1));
                const double wCenterAtKm1 = 0.5 * (wAt(li, lj, lk - 1) + wHere);
                const double dwwdz = (wCenterAtK * wCenterAtK - wCenterAtKm1 * wCenterAtKm1) / dz_;

                const double wEdgeIp = 0.5 * (wHere + wAt(li + 1, lj, lk));
                const double wEdgeIm = 0.5 * (wAt(li - 1, lj, lk) + wHere);
                const double uEdgeKIp = 0.5 * (uAt(li + 1, lj, lk - 1) + uAt(li + 1, lj, lk));
                const double uEdgeKI = 0.5 * (uAt(li, lj, lk - 1) + uAt(li, lj, lk));
                const double dwudx = (wEdgeIp * uEdgeKIp - wEdgeIm * uEdgeKI) / dx_;

                const double wEdgeJp = 0.5 * (wHere + wAt(li, lj + 1, lk));
                const double wEdgeJm = 0.5 * (wAt(li, lj - 1, lk) + wHere);
                const double vEdgeKJp = 0.5 * (vAt(li, lj + 1, lk - 1) + vAt(li, lj + 1, lk));
                const double vEdgeKJ = 0.5 * (vAt(li, lj, lk - 1) + vAt(li, lj, lk));
                const double dwvdy = (wEdgeJp * vEdgeKJp - wEdgeJm * vEdgeKJ) / dy_;

                const double gammaF = viscosity_ + nutAt(li, lj, lk);
                const double gammaB = viscosity_ + nutAt(li, lj, lk - 1);
                const double gammaTransverse = viscosity_ + 0.5 * (nutAt(li, lj, lk - 1) + nutAt(li, lj, lk));

                const double diffusionW =
                    (gammaF * (wAt(li, lj, lk + 1) - wHere) - gammaB * (wHere - wAt(li, lj, lk - 1))) /
                        (dz_ * dz_) +
                    gammaTransverse * (wAt(li + 1, lj, lk) - 2.0 * wHere + wAt(li - 1, lj, lk)) / (dx_ * dx_) +
                    gammaTransverse * (wAt(li, lj + 1, lk) - 2.0 * wHere + wAt(li, lj - 1, lk)) / (dy_ * dy_);

                wStar[indexW(i, j, k)] = wHere + dt * (-(dwwdz + dwudx + dwvdy) + diffusionW);
            }
        }
    }

    // --- k and omega transport: colocated 6-neighbor stencil ---
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

                // --- k transport: destruction beta_star*k*omega*F_DES,
                // Patankar-linearized using the current omega (and the
                // already-computed F_DES) as the frozen factor. F_DES is
                // the one DES-specific change from plain SST -- see the
                // class header.
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
                    const double newK = explicitK / (1.0 + dt * kBetaStar * omega_[idx] * fDes_[idx]);
                    kStar[idx] = std::max(newK, kFloor);
                }

                // --- omega transport: unchanged from plain SST (DES only
                // touches the k-equation's destruction length scale). ---
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
