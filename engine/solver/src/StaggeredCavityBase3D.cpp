#include "aether/solver/ExplicitTimeStep.hpp"
#include "aether/solver/StaggeredCavityBase3D.hpp"

#include "aether/solver/ConvectionLimiter.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace aether::solver {

StaggeredCavityBase3D::StaggeredCavityBase3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                               double lengthY, double lengthZ, double viscosity,
                                               double lidVelocity, ConvectionScheme convection)
    : nx_(nx), ny_(ny), nz_(nz), lengthX_(lengthX), lengthY_(lengthY), lengthZ_(lengthZ),
      dx_(lengthX / static_cast<double>(nx)), dy_(lengthY / static_cast<double>(ny)),
      dz_(lengthZ / static_cast<double>(nz)), viscosity_(viscosity), lidVelocity_(lidVelocity),
      u_((nx + 1) * ny * nz, 0.0), v_(nx * (ny + 1) * nz, 0.0), w_(nx * ny * (nz + 1), 0.0),
      p_(nx * ny * nz, 0.0), convection_(convection) {}

double StaggeredCavityBase3D::schemeTransportValue(double centralValue, double convectingVelocity,
                                                    double near0, double near1, double far0,
                                                    double far1) const {
    if (convection_ == ConvectionScheme::Central) {
        return centralValue;
    }
    const double upwind = convectingVelocity >= 0.0 ? near0 : near1;
    const double downwind = convectingVelocity >= 0.0 ? near1 : near0;
    if (convection_ == ConvectionScheme::FirstOrderUpwind) {
        return upwind;
    }
    const double farUpwind = convectingVelocity >= 0.0 ? far0 : far1;
    const double difference = downwind - upwind;
    if (faceDifferenceIsNegligible(difference, upwind, downwind)) {
        return upwind;
    }
    // Same ratio and blend as LidDrivenCavitySolver2D::schemeFaceValue --
    // see ConvectionLimiter.hpp for why a uniform-grid central estimate
    // makes the classic and gradient-based ratio formulas the same formula.
    const double ratio = (upwind - farUpwind) / difference;
    return upwind + vanLeerLimiter(ratio) * (0.5 * (upwind + downwind) - upwind);
}

void StaggeredCavityBase3D::loadState(std::vector<double> u, std::vector<double> v, std::vector<double> w,
                                       std::vector<double> p, double time) {
    if (u.size() != u_.size() || v.size() != v_.size() || w.size() != w_.size() || p.size() != p_.size()) {
        throw std::invalid_argument("StaggeredCavityBase3D::loadState: field size does not match the grid");
    }
    u_ = std::move(u);
    v_ = std::move(v);
    w_ = std::move(w);
    p_ = std::move(p);
    time_ = time;
}

double StaggeredCavityBase3D::lidVelocityAt(double x, double y) const {
    constexpr double kPi = 3.14159265358979323846;
    const double sx = std::sin(kPi * x / lengthX_);
    const double sy = std::sin(kPi * y / lengthY_);
    return lidVelocity_ * sx * sx * sy * sy;
}

double StaggeredCavityBase3D::wallDistanceAt(std::size_t i, std::size_t j, std::size_t k) const {
    const double x = (static_cast<double>(i) + 0.5) * dx_;
    const double y = (static_cast<double>(j) + 0.5) * dy_;
    const double z = (static_cast<double>(k) + 0.5) * dz_;
    return std::min({x, lengthX_ - x, y, lengthY_ - y, z, lengthZ_ - z});
}

double StaggeredCavityBase3D::nutAt(long long i, long long j, long long k) const {
    if (eddyViscosity_ == nullptr) {
        return 0.0;
    }
    if (i < 0 || i >= static_cast<long long>(nx_) || j < 0 || j >= static_cast<long long>(ny_) || k < 0 ||
        k >= static_cast<long long>(nz_)) {
        return 0.0;
    }
    return (*eddyViscosity_)[indexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j),
                                     static_cast<std::size_t>(k))];
}

double StaggeredCavityBase3D::maxEddyViscosity() const {
    if (eddyViscosity_ == nullptr) {
        return 0.0;
    }
    double maxNut = 0.0;
    for (double n : *eddyViscosity_) {
        maxNut = std::max(maxNut, n);
    }
    return maxNut;
}

double StaggeredCavityBase3D::stableTimeStep() const {
    const double effectiveViscosity = viscosity_ + maxEddyViscosity();
    return explicitStableTimeStep(effectiveViscosity, lidVelocity_, {dx_, dy_, dz_});
}

double StaggeredCavityBase3D::uAt(long long i, long long j, long long k) const {
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

double StaggeredCavityBase3D::vAt(long long i, long long j, long long k) const {
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

double StaggeredCavityBase3D::wAt(long long i, long long j, long long k) const {
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

double StaggeredCavityBase3D::pAt(long long i, long long j, long long k) const {
    const long long ci = std::clamp<long long>(i, 0, static_cast<long long>(nx_) - 1);
    const long long cj = std::clamp<long long>(j, 0, static_cast<long long>(ny_) - 1);
    const long long ck = std::clamp<long long>(k, 0, static_cast<long long>(nz_) - 1);
    return p_[indexP(static_cast<std::size_t>(ci), static_cast<std::size_t>(cj), static_cast<std::size_t>(ck))];
}

double StaggeredCavityBase3D::maxDivergence() const {
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

void StaggeredCavityBase3D::computeMomentumPredictor(std::vector<double>& uStar, std::vector<double>& vStar,
                                                       std::vector<double>& wStar, double dt) const {
    // --- u-momentum: interior x-faces only (i in [1, nx-1]) ---
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 1; i < nx_; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);
                const auto nxLL = static_cast<long long>(nx_);

                const double uHere = uAt(li, lj, lk);
                const double uAtIp1 = uAt(li + 1, lj, lk);
                const double uAtIm1 = uAt(li - 1, lj, lk);
                const double uCenterAtI = 0.5 * (uHere + uAtIp1);
                const double uCenterAtIm1 = 0.5 * (uAtIm1 + uHere);
                // i is u's own staggered direction: uAt has no ghost beyond
                // its two physical boundary faces there (see the class
                // comment), so a stencil point two faces out is clamped to
                // the nearer of those two rather than fetched unchecked.
                const double uAtIp2 = uAt(std::clamp<long long>(li + 2, 0LL, nxLL), lj, lk);
                const double uAtIm2 = uAt(std::clamp<long long>(li - 2, 0LL, nxLL), lj, lk);
                const double transportAtI =
                    schemeTransportValue(uCenterAtI, uCenterAtI, uHere, uAtIp1, uAtIm1, uAtIp2);
                const double transportAtIm1 =
                    schemeTransportValue(uCenterAtIm1, uCenterAtIm1, uAtIm1, uHere, uAtIm2, uAtIp1);
                const double duudx = (uCenterAtI * transportAtI - uCenterAtIm1 * transportAtIm1) / dx_;

                const double uAtJp1 = uAt(li, lj + 1, lk);
                const double uAtJm1 = uAt(li, lj - 1, lk);
                const double uEdgeJp = 0.5 * (uHere + uAtJp1);
                const double uEdgeJm = 0.5 * (uAtJm1 + uHere);
                const double vEdgeIJp = 0.5 * (vAt(li - 1, lj + 1, lk) + vAt(li, lj + 1, lk));
                const double vEdgeIJ = 0.5 * (vAt(li - 1, lj, lk) + vAt(li, lj, lk));
                // j is tangential to u's own face here, so uAt's existing
                // wall mirror already covers any offset safely -- no clamp.
                const double uAtJp2 = uAt(li, lj + 2, lk);
                const double uAtJm2 = uAt(li, lj - 2, lk);
                const double transportJp = schemeTransportValue(uEdgeJp, vEdgeIJp, uHere, uAtJp1, uAtJm1, uAtJp2);
                const double transportJm = schemeTransportValue(uEdgeJm, vEdgeIJ, uAtJm1, uHere, uAtJm2, uAtJp1);
                const double duvdy = (vEdgeIJp * transportJp - vEdgeIJ * transportJm) / dy_;

                const double uAtKp1 = uAt(li, lj, lk + 1);
                const double uAtKm1 = uAt(li, lj, lk - 1);
                const double uEdgeKp = 0.5 * (uHere + uAtKp1);
                const double uEdgeKm = 0.5 * (uAtKm1 + uHere);
                const double wEdgeIKp = 0.5 * (wAt(li - 1, lj, lk + 1) + wAt(li, lj, lk + 1));
                const double wEdgeIK = 0.5 * (wAt(li - 1, lj, lk) + wAt(li, lj, lk));
                const double uAtKp2 = uAt(li, lj, lk + 2);
                const double uAtKm2 = uAt(li, lj, lk - 2);
                const double transportKp = schemeTransportValue(uEdgeKp, wEdgeIKp, uHere, uAtKp1, uAtKm1, uAtKp2);
                const double transportKm = schemeTransportValue(uEdgeKm, wEdgeIK, uAtKm1, uHere, uAtKm2, uAtKp1);
                const double duwdz = (wEdgeIKp * transportKp - wEdgeIK * transportKm) / dz_;

                const double gammaE = viscosity_ + nutAt(li, lj, lk);
                const double gammaW = viscosity_ + nutAt(li - 1, lj, lk);
                const double gammaTransverse = viscosity_ + 0.5 * (nutAt(li - 1, lj, lk) + nutAt(li, lj, lk));

                const double diffusionU =
                    (gammaE * (uAtIp1 - uHere) - gammaW * (uHere - uAtIm1)) / (dx_ * dx_) +
                    gammaTransverse * (uAtJp1 - 2.0 * uHere + uAtJm1) / (dy_ * dy_) +
                    gammaTransverse * (uAtKp1 - 2.0 * uHere + uAtKm1) / (dz_ * dz_);

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
                const auto nyLL = static_cast<long long>(ny_);

                const double vHere = vAt(li, lj, lk);
                const double vAtJp1 = vAt(li, lj + 1, lk);
                const double vAtJm1 = vAt(li, lj - 1, lk);
                const double vCenterAtJ = 0.5 * (vHere + vAtJp1);
                const double vCenterAtJm1 = 0.5 * (vAtJm1 + vHere);
                // j is v's own staggered direction -- same clamp as u's own
                // direction above, and for the same reason.
                const double vAtJp2 = vAt(li, std::clamp<long long>(lj + 2, 0LL, nyLL), lk);
                const double vAtJm2 = vAt(li, std::clamp<long long>(lj - 2, 0LL, nyLL), lk);
                const double transportAtJ =
                    schemeTransportValue(vCenterAtJ, vCenterAtJ, vHere, vAtJp1, vAtJm1, vAtJp2);
                const double transportAtJm1 =
                    schemeTransportValue(vCenterAtJm1, vCenterAtJm1, vAtJm1, vHere, vAtJm2, vAtJp1);
                const double dvvdy = (vCenterAtJ * transportAtJ - vCenterAtJm1 * transportAtJm1) / dy_;

                const double vAtKp1 = vAt(li, lj, lk + 1);
                const double vAtKm1 = vAt(li, lj, lk - 1);
                const double vEdgeKp = 0.5 * (vHere + vAtKp1);
                const double vEdgeKm = 0.5 * (vAtKm1 + vHere);
                const double wEdgeJKp = 0.5 * (wAt(li, lj - 1, lk + 1) + wAt(li, lj, lk + 1));
                const double wEdgeJK = 0.5 * (wAt(li, lj - 1, lk) + wAt(li, lj, lk));
                const double vAtKp2 = vAt(li, lj, lk + 2);
                const double vAtKm2 = vAt(li, lj, lk - 2);
                const double transportKp = schemeTransportValue(vEdgeKp, wEdgeJKp, vHere, vAtKp1, vAtKm1, vAtKp2);
                const double transportKm = schemeTransportValue(vEdgeKm, wEdgeJK, vAtKm1, vHere, vAtKm2, vAtKp1);
                const double dvwdz = (wEdgeJKp * transportKp - wEdgeJK * transportKm) / dz_;

                const double vAtIp1 = vAt(li + 1, lj, lk);
                const double vAtIm1 = vAt(li - 1, lj, lk);
                const double vEdgeIp = 0.5 * (vHere + vAtIp1);
                const double vEdgeIm = 0.5 * (vAtIm1 + vHere);
                const double uEdgeJIp = 0.5 * (uAt(li + 1, lj - 1, lk) + uAt(li + 1, lj, lk));
                const double uEdgeJI = 0.5 * (uAt(li, lj - 1, lk) + uAt(li, lj, lk));
                const double vAtIp2 = vAt(li + 2, lj, lk);
                const double vAtIm2 = vAt(li - 2, lj, lk);
                const double transportIp = schemeTransportValue(vEdgeIp, uEdgeJIp, vHere, vAtIp1, vAtIm1, vAtIp2);
                const double transportIm = schemeTransportValue(vEdgeIm, uEdgeJI, vAtIm1, vHere, vAtIm2, vAtIp1);
                const double dvudx = (vEdgeIp * transportIp - vEdgeIm * transportIm) / dx_;

                const double gammaN = viscosity_ + nutAt(li, lj, lk);
                const double gammaS = viscosity_ + nutAt(li, lj - 1, lk);
                const double gammaTransverse = viscosity_ + 0.5 * (nutAt(li, lj - 1, lk) + nutAt(li, lj, lk));

                const double diffusionV =
                    (gammaN * (vAtJp1 - vHere) - gammaS * (vHere - vAtJm1)) / (dy_ * dy_) +
                    gammaTransverse * (vAtIp1 - 2.0 * vHere + vAtIm1) / (dx_ * dx_) +
                    gammaTransverse * (vAtKp1 - 2.0 * vHere + vAtKm1) / (dz_ * dz_);

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
                const auto nzLL = static_cast<long long>(nz_);

                const double wHere = wAt(li, lj, lk);
                const double wAtKp1 = wAt(li, lj, lk + 1);
                const double wAtKm1 = wAt(li, lj, lk - 1);
                const double wCenterAtK = 0.5 * (wHere + wAtKp1);
                const double wCenterAtKm1 = 0.5 * (wAtKm1 + wHere);
                // k is w's own staggered direction -- same clamp as u's and
                // v's own directions above, and for the same reason.
                const double wAtKp2 = wAt(li, lj, std::clamp<long long>(lk + 2, 0LL, nzLL));
                const double wAtKm2 = wAt(li, lj, std::clamp<long long>(lk - 2, 0LL, nzLL));
                const double transportAtK =
                    schemeTransportValue(wCenterAtK, wCenterAtK, wHere, wAtKp1, wAtKm1, wAtKp2);
                const double transportAtKm1 =
                    schemeTransportValue(wCenterAtKm1, wCenterAtKm1, wAtKm1, wHere, wAtKm2, wAtKp1);
                const double dwwdz = (wCenterAtK * transportAtK - wCenterAtKm1 * transportAtKm1) / dz_;

                const double wAtIp1 = wAt(li + 1, lj, lk);
                const double wAtIm1 = wAt(li - 1, lj, lk);
                const double wEdgeIp = 0.5 * (wHere + wAtIp1);
                const double wEdgeIm = 0.5 * (wAtIm1 + wHere);
                const double uEdgeKIp = 0.5 * (uAt(li + 1, lj, lk - 1) + uAt(li + 1, lj, lk));
                const double uEdgeKI = 0.5 * (uAt(li, lj, lk - 1) + uAt(li, lj, lk));
                const double wAtIp2 = wAt(li + 2, lj, lk);
                const double wAtIm2 = wAt(li - 2, lj, lk);
                const double transportIp = schemeTransportValue(wEdgeIp, uEdgeKIp, wHere, wAtIp1, wAtIm1, wAtIp2);
                const double transportIm = schemeTransportValue(wEdgeIm, uEdgeKI, wAtIm1, wHere, wAtIm2, wAtIp1);
                const double dwudx = (wEdgeIp * transportIp - wEdgeIm * transportIm) / dx_;

                const double wAtJp1 = wAt(li, lj + 1, lk);
                const double wAtJm1 = wAt(li, lj - 1, lk);
                const double wEdgeJp = 0.5 * (wHere + wAtJp1);
                const double wEdgeJm = 0.5 * (wAtJm1 + wHere);
                const double vEdgeKJp = 0.5 * (vAt(li, lj + 1, lk - 1) + vAt(li, lj + 1, lk));
                const double vEdgeKJ = 0.5 * (vAt(li, lj, lk - 1) + vAt(li, lj, lk));
                const double wAtJp2 = wAt(li, lj + 2, lk);
                const double wAtJm2 = wAt(li, lj - 2, lk);
                const double transportJp = schemeTransportValue(wEdgeJp, vEdgeKJp, wHere, wAtJp1, wAtJm1, wAtJp2);
                const double transportJm = schemeTransportValue(wEdgeJm, vEdgeKJ, wAtJm1, wHere, wAtJm2, wAtJp1);
                const double dwvdy = (vEdgeKJp * transportJp - vEdgeKJ * transportJm) / dy_;

                const double gammaF = viscosity_ + nutAt(li, lj, lk);
                const double gammaB = viscosity_ + nutAt(li, lj, lk - 1);
                const double gammaTransverse = viscosity_ + 0.5 * (nutAt(li, lj, lk - 1) + nutAt(li, lj, lk));

                const double diffusionW =
                    (gammaF * (wAtKp1 - wHere) - gammaB * (wHere - wAtKm1)) / (dz_ * dz_) +
                    gammaTransverse * (wAtIp1 - 2.0 * wHere + wAtIm1) / (dx_ * dx_) +
                    gammaTransverse * (wAtJp1 - 2.0 * wHere + wAtJm1) / (dy_ * dy_);

                wStar[indexW(i, j, k)] = wHere + dt * (-(dwwdz + dwudx + dwvdy) + diffusionW);
            }
        }
    }
}

std::vector<double> StaggeredCavityBase3D::applyLaplacian(const std::vector<double>& x) const {
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

double StaggeredCavityBase3D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void StaggeredCavityBase3D::projectToDivergenceFree(std::vector<double>& uStar, std::vector<double>& vStar,
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

} // namespace aether::solver
