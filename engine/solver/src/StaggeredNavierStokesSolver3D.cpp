#include "aether/solver/ExplicitTimeStep.hpp"
#include "aether/solver/StaggeredNavierStokesSolver3D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

StaggeredNavierStokesSolver3D::StaggeredNavierStokesSolver3D(std::size_t nx, std::size_t ny,
                                                               std::size_t nz, double lengthX,
                                                               double lengthY, double lengthZ,
                                                               double viscosity)
    : nx_(nx), ny_(ny), nz_(nz), dx_(lengthX / static_cast<double>(nx)),
      dy_(lengthY / static_cast<double>(ny)), dz_(lengthZ / static_cast<double>(nz)),
      viscosity_(viscosity), u_(nx * ny * nz, 0.0), v_(nx * ny * nz, 0.0), w_(nx * ny * nz, 0.0),
      p_(nx * ny * nz, 0.0) {}

void StaggeredNavierStokesSolver3D::setVelocity(std::size_t i, std::size_t j, std::size_t k, double u,
                                                 double v, double w) {
    const std::size_t idx = index(i, j, k);
    u_[idx] = u;
    v_[idx] = v;
    w_[idx] = w;
}

std::size_t StaggeredNavierStokesSolver3D::wrap(long long i, std::size_t n) const {
    long long m = i % static_cast<long long>(n);
    if (m < 0) {
        m += static_cast<long long>(n);
    }
    return static_cast<std::size_t>(m);
}

double StaggeredNavierStokesSolver3D::stableTimeStep(double velocityScale) const {
    return explicitStableTimeStep(viscosity_, velocityScale, {dx_, dy_, dz_});
}

double StaggeredNavierStokesSolver3D::u(std::size_t i, std::size_t j, std::size_t k) const {
    return u_[index(i, j, k)];
}
double StaggeredNavierStokesSolver3D::v(std::size_t i, std::size_t j, std::size_t k) const {
    return v_[index(i, j, k)];
}
double StaggeredNavierStokesSolver3D::w(std::size_t i, std::size_t j, std::size_t k) const {
    return w_[index(i, j, k)];
}
double StaggeredNavierStokesSolver3D::pressure(std::size_t i, std::size_t j, std::size_t k) const {
    return p_[index(i, j, k)];
}

double StaggeredNavierStokesSolver3D::maxDivergence() const {
    double maxDiv = 0.0;
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = index(i, j, k);
                const std::size_t ip = index(wrap(static_cast<long long>(i) + 1, nx_), j, k);
                const std::size_t jp = index(i, wrap(static_cast<long long>(j) + 1, ny_), k);
                const std::size_t kp = index(i, j, wrap(static_cast<long long>(k) + 1, nz_));
                const double div = (u_[ip] - u_[idx]) / dx_ + (v_[jp] - v_[idx]) / dy_ +
                                    (w_[kp] - w_[idx]) / dz_;
                maxDiv = std::max(maxDiv, std::fabs(div));
            }
        }
    }
    return maxDiv;
}

std::vector<double>
StaggeredNavierStokesSolver3D::applyPeriodicLaplacian(const std::vector<double>& x) const {
    const double ax = 1.0 / (dx_ * dx_);
    const double ay = 1.0 / (dy_ * dy_);
    const double az = 1.0 / (dz_ * dz_);
    const double weightTotal = 2.0 * ax + 2.0 * ay + 2.0 * az;

    std::vector<double> result(x.size());
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = index(i, j, k);
                if (idx == 0) {
                    result[idx] = x[idx]; // pinned reference cell: removes the null space
                    continue;
                }
                const std::size_t ip = index(wrap(static_cast<long long>(i) + 1, nx_), j, k);
                const std::size_t im = index(wrap(static_cast<long long>(i) - 1, nx_), j, k);
                const std::size_t jp = index(i, wrap(static_cast<long long>(j) + 1, ny_), k);
                const std::size_t jm = index(i, wrap(static_cast<long long>(j) - 1, ny_), k);
                const std::size_t kp = index(i, j, wrap(static_cast<long long>(k) + 1, nz_));
                const std::size_t km = index(i, j, wrap(static_cast<long long>(k) - 1, nz_));
                const double weightedSum =
                    ax * (x[ip] + x[im]) + ay * (x[jp] + x[jm]) + az * (x[kp] + x[km]);
                result[idx] = weightTotal * x[idx] - weightedSum;
            }
        }
    }
    return result;
}

double StaggeredNavierStokesSolver3D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void StaggeredNavierStokesSolver3D::projectToDivergenceFree(std::vector<double>& uStar,
                                                             std::vector<double>& vStar,
                                                             std::vector<double>& wStar, double dt) {
    const std::size_t n = uStar.size();

    std::vector<double> rhs(n);
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = index(i, j, k);
                if (idx == 0) {
                    rhs[idx] = 0.0;
                    continue;
                }
                const std::size_t ip = index(wrap(static_cast<long long>(i) + 1, nx_), j, k);
                const std::size_t jp = index(i, wrap(static_cast<long long>(j) + 1, ny_), k);
                const std::size_t kp = index(i, j, wrap(static_cast<long long>(k) + 1, nz_));
                const double divergence = (uStar[ip] - uStar[idx]) / dx_ +
                                           (vStar[jp] - vStar[idx]) / dy_ +
                                           (wStar[kp] - wStar[idx]) / dz_;
                // applyPeriodicLaplacian() computes -nabla^2(p), so the
                // right-hand side of nabla^2(p) = div(uStar)/dt must be
                // negated (same convention as TaylorGreenVortexSolver2D).
                rhs[idx] = -divergence / dt;
            }
        }
    }

    // Conjugate Gradient on the periodic pressure-Poisson system, same
    // pattern as TaylorGreenVortexSolver2D::projectToDivergenceFree().
    std::vector<double> residual(n);
    {
        const std::vector<double> operatorAppliedToP = applyPeriodicLaplacian(p_);
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
        const std::vector<double> operatorAppliedToDirection = applyPeriodicLaplacian(direction);
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

    // Correct the predicted velocity with the solved pressure's gradient,
    // each component at its own staggered face location (adjacent cell
    // pair, h-wide -- not the collocated solvers' 2h-wide central
    // difference).
    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const std::size_t idx = index(i, j, k);
                const std::size_t im = index(wrap(static_cast<long long>(i) - 1, nx_), j, k);
                const std::size_t jm = index(i, wrap(static_cast<long long>(j) - 1, ny_), k);
                const std::size_t km = index(i, j, wrap(static_cast<long long>(k) - 1, nz_));
                const double dpdx = (p_[idx] - p_[im]) / dx_;
                const double dpdy = (p_[idx] - p_[jm]) / dy_;
                const double dpdz = (p_[idx] - p_[km]) / dz_;
                u_[idx] = uStar[idx] - dt * dpdx;
                v_[idx] = vStar[idx] - dt * dpdy;
                w_[idx] = wStar[idx] - dt * dpdz;
            }
        }
    }
}

void StaggeredNavierStokesSolver3D::step(double dt) {
    const std::size_t n = u_.size();
    std::vector<double> uStar(n);
    std::vector<double> vStar(n);
    std::vector<double> wStar(n);

    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                // Scalar wrapped coordinates -- for building *new* (i,j,k)
                // triples passed back into index(), e.g. index(iM, jP, k).
                const std::size_t iP = wrap(static_cast<long long>(i) + 1, nx_);
                const std::size_t iM = wrap(static_cast<long long>(i) - 1, nx_);
                const std::size_t jP = wrap(static_cast<long long>(j) + 1, ny_);
                const std::size_t jM = wrap(static_cast<long long>(j) - 1, ny_);
                const std::size_t kP = wrap(static_cast<long long>(k) + 1, nz_);
                const std::size_t kM = wrap(static_cast<long long>(k) - 1, nz_);

                // Flattened neighbor indices -- for direct same-field,
                // same-axis lookups only (u_[ip], v_[jm], etc.). These must
                // never be reused as a scalar coordinate argument to
                // index(): they are already a full flattened index, not a
                // single-axis position.
                const std::size_t idx = index(i, j, k);
                const std::size_t ip = index(iP, j, k);
                const std::size_t im = index(iM, j, k);
                const std::size_t jp = index(i, jP, k);
                const std::size_t jm = index(i, jM, k);
                const std::size_t kp = index(i, j, kP);
                const std::size_t km = index(i, j, kM);

                // --- u-momentum, at this u-face (i*dx, (j+0.5)dy, (k+0.5)dz) ---
                {
                    const double uCenterAtI = 0.5 * (u_[idx] + u_[ip]);
                    const double uCenterAtIm1 = 0.5 * (u_[im] + u_[idx]);
                    const double duudx = (uCenterAtI * uCenterAtI - uCenterAtIm1 * uCenterAtIm1) / dx_;

                    const double uEdgeJp = 0.5 * (u_[idx] + u_[jp]);
                    const double uEdgeJm = 0.5 * (u_[jm] + u_[idx]);
                    const double vEdgeIJp = 0.5 * (v_[index(iM, jP, k)] + v_[index(i, jP, k)]);
                    const double vEdgeIJ = 0.5 * (v_[index(iM, j, k)] + v_[idx]);
                    const double duvdy = (uEdgeJp * vEdgeIJp - uEdgeJm * vEdgeIJ) / dy_;

                    const double uEdgeKp = 0.5 * (u_[idx] + u_[kp]);
                    const double uEdgeKm = 0.5 * (u_[km] + u_[idx]);
                    const double wEdgeIKp = 0.5 * (w_[index(iM, j, kP)] + w_[index(i, j, kP)]);
                    const double wEdgeIK = 0.5 * (w_[index(iM, j, k)] + w_[idx]);
                    const double duwdz = (uEdgeKp * wEdgeIKp - uEdgeKm * wEdgeIK) / dz_;

                    const double laplacianU = (u_[ip] - 2.0 * u_[idx] + u_[im]) / (dx_ * dx_) +
                                               (u_[jp] - 2.0 * u_[idx] + u_[jm]) / (dy_ * dy_) +
                                               (u_[kp] - 2.0 * u_[idx] + u_[km]) / (dz_ * dz_);

                    uStar[idx] = u_[idx] + dt * (-(duudx + duvdy + duwdz) + viscosity_ * laplacianU);
                }

                // --- v-momentum, at this v-face ((i+0.5)dx, j*dy, (k+0.5)dz) ---
                {
                    const double vCenterAtJ = 0.5 * (v_[idx] + v_[jp]);
                    const double vCenterAtJm1 = 0.5 * (v_[jm] + v_[idx]);
                    const double dvvdy = (vCenterAtJ * vCenterAtJ - vCenterAtJm1 * vCenterAtJm1) / dy_;

                    const double vEdgeKp = 0.5 * (v_[idx] + v_[kp]);
                    const double vEdgeKm = 0.5 * (v_[km] + v_[idx]);
                    const double wEdgeJKp = 0.5 * (w_[index(i, jM, kP)] + w_[index(i, j, kP)]);
                    const double wEdgeJK = 0.5 * (w_[index(i, jM, k)] + w_[idx]);
                    const double dvwdz = (vEdgeKp * wEdgeJKp - vEdgeKm * wEdgeJK) / dz_;

                    const double vEdgeIp = 0.5 * (v_[idx] + v_[ip]);
                    const double vEdgeIm = 0.5 * (v_[im] + v_[idx]);
                    const double uEdgeJIp = 0.5 * (u_[index(iP, jM, k)] + u_[index(iP, j, k)]);
                    const double uEdgeJI = 0.5 * (u_[index(i, jM, k)] + u_[idx]);
                    const double dvudx = (vEdgeIp * uEdgeJIp - vEdgeIm * uEdgeJI) / dx_;

                    const double laplacianV = (v_[ip] - 2.0 * v_[idx] + v_[im]) / (dx_ * dx_) +
                                               (v_[jp] - 2.0 * v_[idx] + v_[jm]) / (dy_ * dy_) +
                                               (v_[kp] - 2.0 * v_[idx] + v_[km]) / (dz_ * dz_);

                    vStar[idx] = v_[idx] + dt * (-(dvvdy + dvwdz + dvudx) + viscosity_ * laplacianV);
                }

                // --- w-momentum, at this w-face ((i+0.5)dx, (j+0.5)dy, k*dz) ---
                {
                    const double wCenterAtK = 0.5 * (w_[idx] + w_[kp]);
                    const double wCenterAtKm1 = 0.5 * (w_[km] + w_[idx]);
                    const double dwwdz = (wCenterAtK * wCenterAtK - wCenterAtKm1 * wCenterAtKm1) / dz_;

                    const double wEdgeIp = 0.5 * (w_[idx] + w_[ip]);
                    const double wEdgeIm = 0.5 * (w_[im] + w_[idx]);
                    const double uEdgeKIp = 0.5 * (u_[index(iP, j, kM)] + u_[index(iP, j, k)]);
                    const double uEdgeKI = 0.5 * (u_[index(i, j, kM)] + u_[idx]);
                    const double dwudx = (wEdgeIp * uEdgeKIp - wEdgeIm * uEdgeKI) / dx_;

                    const double wEdgeJp = 0.5 * (w_[idx] + w_[jp]);
                    const double wEdgeJm = 0.5 * (w_[jm] + w_[idx]);
                    const double vEdgeKJp = 0.5 * (v_[index(i, jP, kM)] + v_[index(i, jP, k)]);
                    const double vEdgeKJ = 0.5 * (v_[index(i, j, kM)] + v_[idx]);
                    const double dwvdy = (wEdgeJp * vEdgeKJp - wEdgeJm * vEdgeKJ) / dy_;

                    const double laplacianW = (w_[ip] - 2.0 * w_[idx] + w_[im]) / (dx_ * dx_) +
                                               (w_[jp] - 2.0 * w_[idx] + w_[jm]) / (dy_ * dy_) +
                                               (w_[kp] - 2.0 * w_[idx] + w_[km]) / (dz_ * dz_);

                    wStar[idx] = w_[idx] + dt * (-(dwwdz + dwudx + dwvdy) + viscosity_ * laplacianW);
                }
            }
        }
    }

    projectToDivergenceFree(uStar, vStar, wStar, dt);
    time_ += dt;
}

} // namespace aether::solver
