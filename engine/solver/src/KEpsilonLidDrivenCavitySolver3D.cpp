#include "aether/solver/KEpsilonLidDrivenCavitySolver3D.hpp"

#include "aether/solver/MixingLengthLidDrivenCavitySolver3D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

namespace {
constexpr double kCmu = 0.09;
constexpr double kCeps1 = 1.44;
constexpr double kCeps2 = 1.92;
constexpr double kSigmaK = 1.0;
constexpr double kSigmaEps = 1.3;
constexpr double kFloor = 1e-10;
} // namespace

KEpsilonLidDrivenCavitySolver3D::KEpsilonLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny,
                                                                   std::size_t nz, double lengthX,
                                                                   double lengthY, double lengthZ,
                                                                   double viscosity, double lidVelocity)
    : StaggeredCavityBase3D(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity),
      k_(nx * ny * nz, 0.0), epsilon_(nx * ny * nz, 0.0), nut_(nx * ny * nz, 0.0),
      production_(nx * ny * nz, 0.0) {
    // Warm-start velocity from the simpler mixing-length closure -- same
    // cold-start deadlock reason the 2D class warm-starts from
    // MixingLengthLidDrivenCavitySolver2D (u=0 gives zero production, so
    // k/epsilon would collapse before the velocity field ever develops).
    // The base's momentum predictor and dt limit read nu_t through this.
    setEddyViscosityField(&nut_);

    MixingLengthLidDrivenCavitySolver3D primer(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity);
    const double primerDt = primer.stableTimeStep();
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

    // Small, uniform initial k/epsilon -- same self-consistent starting
    // guess as the 2D class (eps0 from k0 via eps ~ Cmu^0.75 * k^1.5 / L).
    const double k0 = 0.01 * lidVelocity * lidVelocity;
    const double lengthScale = 0.1 * std::min({lengthX_, lengthY_, lengthZ_});
    const double eps0 = std::pow(kCmu, 0.75) * std::pow(k0, 1.5) / lengthScale;
    std::fill(k_.begin(), k_.end(), k0);
    std::fill(epsilon_.begin(), epsilon_.end(), eps0);
}

double KEpsilonLidDrivenCavitySolver3D::k(std::size_t i, std::size_t j, std::size_t k) const {
    return k_[indexP(i, j, k)];
}

double KEpsilonLidDrivenCavitySolver3D::epsilon(std::size_t i, std::size_t j, std::size_t k) const {
    return epsilon_[indexP(i, j, k)];
}

double KEpsilonLidDrivenCavitySolver3D::eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const {
    return nut_[indexP(i, j, k)];
}

double KEpsilonLidDrivenCavitySolver3D::kAt(long long i, long long j, long long k) const {
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

double KEpsilonLidDrivenCavitySolver3D::epsilonGhostAt(std::size_t i, std::size_t j, std::size_t k, int di,
                                                         int dj, int dk) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    const long long nk = static_cast<long long>(k) + dk;
    const bool crossesWall = ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 ||
                              nj >= static_cast<long long>(ny_) || nk < 0 || nk >= static_cast<long long>(nz_);
    if (!crossesWall) {
        return epsilon_[indexP(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj),
                                static_cast<std::size_t>(nk))];
    }
    const double halfSpacing = (di != 0 ? dx_ : (dj != 0 ? dy_ : dz_)) / 2.0;
    const double kNearWall = std::max(k_[indexP(i, j, k)], 0.0);
    const double epsWall = 2.0 * viscosity_ * kNearWall / (halfSpacing * halfSpacing);
    return 2.0 * epsWall - epsilon_[indexP(i, j, k)];
}

void KEpsilonLidDrivenCavitySolver3D::updateEddyViscosityAndProduction() {
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = indexP(i, j, k);
                const double kSafe = std::max(k_[idx], kFloor);
                const double epsSafe = std::max(epsilon_[idx], kFloor);
                nut_[idx] = kCmu * kSafe * kSafe / epsSafe;

                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

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

                production_[idx] = nut_[idx] * strainSquared;
            }
        }
    }
}

void KEpsilonLidDrivenCavitySolver3D::step(double dt) {
    updateEddyViscosityAndProduction();

    std::vector<double> uStar = u_;
    std::vector<double> vStar = v_;
    std::vector<double> wStar = w_;
    std::vector<double> kStar(k_.size());
    std::vector<double> epsStar(epsilon_.size());

    computeMomentumPredictor(uStar, vStar, wStar, dt);

    // --- k and epsilon transport: plain colocated 6-neighbor stencil,
    // both fields living at cell centers like pressure and nu_t. ---
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

                // --- k transport ---
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

                    const double gEk = viscosity_ + 0.5 * (nutHere + nutAt(li + 1, lj, lk)) / kSigmaK;
                    const double gWk = viscosity_ + 0.5 * (nutHere + nutAt(li - 1, lj, lk)) / kSigmaK;
                    const double gNk = viscosity_ + 0.5 * (nutHere + nutAt(li, lj + 1, lk)) / kSigmaK;
                    const double gSk = viscosity_ + 0.5 * (nutHere + nutAt(li, lj - 1, lk)) / kSigmaK;
                    const double gFk = viscosity_ + 0.5 * (nutHere + nutAt(li, lj, lk + 1)) / kSigmaK;
                    const double gBk = viscosity_ + 0.5 * (nutHere + nutAt(li, lj, lk - 1)) / kSigmaK;

                    const double diffusionK = (gEk * (kE - kHere) - gWk * (kHere - kW)) / (dx_ * dx_) +
                                               (gNk * (kN - kHere) - gSk * (kHere - kS)) / (dy_ * dy_) +
                                               (gFk * (kF - kHere) - gBk * (kHere - kB)) / (dz_ * dz_);

                    const double newK =
                        kHere + dt * (-advectionK + production_[idx] - epsilon_[idx] + diffusionK);
                    kStar[idx] = std::max(newK, kFloor);
                }

                // --- epsilon transport (Patankar-linearized destruction,
                // same frozen-factor implicit form as the 2D class). ---
                {
                    const double epsE = epsilonGhostAt(i, j, k, 1, 0, 0);
                    const double epsW = epsilonGhostAt(i, j, k, -1, 0, 0);
                    const double epsN = epsilonGhostAt(i, j, k, 0, 1, 0);
                    const double epsS = epsilonGhostAt(i, j, k, 0, -1, 0);
                    const double epsF = epsilonGhostAt(i, j, k, 0, 0, 1);
                    const double epsB = epsilonGhostAt(i, j, k, 0, 0, -1);
                    const double epsHere = epsilon_[idx];

                    const double depsdx = (epsE - epsW) / (2.0 * dx_);
                    const double depsdy = (epsN - epsS) / (2.0 * dy_);
                    const double depsdz = (epsF - epsB) / (2.0 * dz_);
                    const double advectionEps = uCenter * depsdx + vCenter * depsdy + wCenter * depsdz;

                    const double gEe = viscosity_ + 0.5 * (nutHere + nutAt(li + 1, lj, lk)) / kSigmaEps;
                    const double gWe = viscosity_ + 0.5 * (nutHere + nutAt(li - 1, lj, lk)) / kSigmaEps;
                    const double gNe = viscosity_ + 0.5 * (nutHere + nutAt(li, lj + 1, lk)) / kSigmaEps;
                    const double gSe = viscosity_ + 0.5 * (nutHere + nutAt(li, lj - 1, lk)) / kSigmaEps;
                    const double gFe = viscosity_ + 0.5 * (nutHere + nutAt(li, lj, lk + 1)) / kSigmaEps;
                    const double gBe = viscosity_ + 0.5 * (nutHere + nutAt(li, lj, lk - 1)) / kSigmaEps;

                    const double diffusionEps =
                        (gEe * (epsE - epsHere) - gWe * (epsHere - epsW)) / (dx_ * dx_) +
                        (gNe * (epsN - epsHere) - gSe * (epsHere - epsS)) / (dy_ * dy_) +
                        (gFe * (epsF - epsHere) - gBe * (epsHere - epsB)) / (dz_ * dz_);

                    const double kSafe = std::max(k_[idx], kFloor);
                    const double productionEps = kCeps1 * (epsHere / kSafe) * production_[idx];
                    const double destructionRate = kCeps2 * epsHere / kSafe; // frozen factor

                    const double explicitPart =
                        epsHere + dt * (-advectionEps + productionEps + diffusionEps);
                    const double newEps = explicitPart / (1.0 + dt * destructionRate);
                    epsStar[idx] = std::max(newEps, kFloor);
                }
            }
        }
    }

    k_ = std::move(kStar);
    epsilon_ = std::move(epsStar);

    projectToDivergenceFree(uStar, vStar, wStar, dt);
    time_ += dt;
}

} // namespace aether::solver
