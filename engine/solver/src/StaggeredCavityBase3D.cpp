#include "aether/solver/StaggeredCavityBase3D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

StaggeredCavityBase3D::StaggeredCavityBase3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                               double lengthY, double lengthZ, double viscosity,
                                               double lidVelocity)
    : nx_(nx), ny_(ny), nz_(nz), lengthX_(lengthX), lengthY_(lengthY), lengthZ_(lengthZ),
      dx_(lengthX / static_cast<double>(nx)), dy_(lengthY / static_cast<double>(ny)),
      dz_(lengthZ / static_cast<double>(nz)), viscosity_(viscosity), lidVelocity_(lidVelocity),
      u_((nx + 1) * ny * nz, 0.0), v_(nx * (ny + 1) * nz, 0.0), w_(nx * ny * (nz + 1), 0.0),
      p_(nx * ny * nz, 0.0) {}

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
    const double diffusiveLimit = 1.0 / (2.0 * effectiveViscosity *
                                          (1.0 / (dx_ * dx_) + 1.0 / (dy_ * dy_) + 1.0 / (dz_ * dz_)));
    const double convectiveLimit = std::min({dx_, dy_, dz_}) / std::max(lidVelocity_, 1e-12);
    return std::min(diffusiveLimit, convectiveLimit);
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
