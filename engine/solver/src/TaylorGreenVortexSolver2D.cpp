#include "aether/solver/ExplicitTimeStep.hpp"
#include "aether/solver/TaylorGreenVortexSolver2D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

TaylorGreenVortexSolver2D::TaylorGreenVortexSolver2D(std::size_t nx, std::size_t ny, double lengthX,
                                                       double lengthY, double viscosity)
    : nx_(nx), ny_(ny), dx_(lengthX / static_cast<double>(nx)), dy_(lengthY / static_cast<double>(ny)),
      viscosity_(viscosity), u_(nx * ny, 0.0), v_(nx * ny, 0.0), p_(nx * ny, 0.0) {}

void TaylorGreenVortexSolver2D::setVelocity(std::size_t i, std::size_t j, double u, double v) {
    const std::size_t idx = index(i, j);
    u_[idx] = u;
    v_[idx] = v;
}

std::size_t TaylorGreenVortexSolver2D::wrap(long long i, std::size_t n) const {
    long long m = i % static_cast<long long>(n);
    if (m < 0) {
        m += static_cast<long long>(n);
    }
    return static_cast<std::size_t>(m);
}

double TaylorGreenVortexSolver2D::stableTimeStep(double velocityScale) const {
    return explicitStableTimeStep(viscosity_, velocityScale, {dx_, dy_});
}

double TaylorGreenVortexSolver2D::u(std::size_t i, std::size_t j) const { return u_[index(i, j)]; }
double TaylorGreenVortexSolver2D::v(std::size_t i, std::size_t j) const { return v_[index(i, j)]; }
double TaylorGreenVortexSolver2D::pressure(std::size_t i, std::size_t j) const { return p_[index(i, j)]; }

double TaylorGreenVortexSolver2D::maxDivergence() const {
    double maxDiv = 0.0;
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t ip = index(wrap(static_cast<long long>(i) + 1, nx_), j);
            const std::size_t im = index(wrap(static_cast<long long>(i) - 1, nx_), j);
            const std::size_t jp = index(i, wrap(static_cast<long long>(j) + 1, ny_));
            const std::size_t jm = index(i, wrap(static_cast<long long>(j) - 1, ny_));
            const double div = (u_[ip] - u_[im]) / (2.0 * dx_) + (v_[jp] - v_[jm]) / (2.0 * dy_);
            maxDiv = std::max(maxDiv, std::fabs(div));
        }
    }
    return maxDiv;
}

std::vector<double> TaylorGreenVortexSolver2D::applyPeriodicLaplacian(const std::vector<double>& x) const {
    const double ax = 1.0 / (dx_ * dx_);
    const double ay = 1.0 / (dy_ * dy_);
    const double weightTotal = 2.0 * ax + 2.0 * ay;

    std::vector<double> result(x.size());
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            if (idx == 0) {
                result[idx] = x[idx]; // pinned reference cell: removes the null space
                continue;
            }
            const std::size_t ip = index(wrap(static_cast<long long>(i) + 1, nx_), j);
            const std::size_t im = index(wrap(static_cast<long long>(i) - 1, nx_), j);
            const std::size_t jp = index(i, wrap(static_cast<long long>(j) + 1, ny_));
            const std::size_t jm = index(i, wrap(static_cast<long long>(j) - 1, ny_));
            const double weightedSum = ax * (x[ip] + x[im]) + ay * (x[jp] + x[jm]);
            result[idx] = weightTotal * x[idx] - weightedSum;
        }
    }
    return result;
}

double TaylorGreenVortexSolver2D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void TaylorGreenVortexSolver2D::projectToDivergenceFree(std::vector<double>& uStar,
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
            const std::size_t ip = index(wrap(static_cast<long long>(i) + 1, nx_), j);
            const std::size_t im = index(wrap(static_cast<long long>(i) - 1, nx_), j);
            const std::size_t jp = index(i, wrap(static_cast<long long>(j) + 1, ny_));
            const std::size_t jm = index(i, wrap(static_cast<long long>(j) - 1, ny_));
            const double divergence =
                (uStar[ip] - uStar[im]) / (2.0 * dx_) + (vStar[jp] - vStar[jm]) / (2.0 * dy_);
            // applyPeriodicLaplacian() computes -nabla^2(p) (matching
            // DiffusionProblem's convention), so to solve the intended
            // nabla^2(p) = div(uStar)/dt, the right-hand side here must be
            // negated.
            rhs[idx] = -divergence / dt;
        }
    }

    // Conjugate Gradient on the periodic pressure-Poisson system, reusing
    // the same pattern as SteadyDiffusionSolver::solveConjugateGradient().
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

    // Correct the predicted velocity with the solved pressure's gradient.
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            const std::size_t ip = index(wrap(static_cast<long long>(i) + 1, nx_), j);
            const std::size_t im = index(wrap(static_cast<long long>(i) - 1, nx_), j);
            const std::size_t jp = index(i, wrap(static_cast<long long>(j) + 1, ny_));
            const std::size_t jm = index(i, wrap(static_cast<long long>(j) - 1, ny_));
            const double dpdx = (p_[ip] - p_[im]) / (2.0 * dx_);
            const double dpdy = (p_[jp] - p_[jm]) / (2.0 * dy_);
            u_[idx] = uStar[idx] - dt * dpdx;
            v_[idx] = vStar[idx] - dt * dpdy;
        }
    }
}

void TaylorGreenVortexSolver2D::step(double dt) {
    std::vector<double> uStar(u_.size());
    std::vector<double> vStar(v_.size());

    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            const std::size_t ip = index(wrap(static_cast<long long>(i) + 1, nx_), j);
            const std::size_t im = index(wrap(static_cast<long long>(i) - 1, nx_), j);
            const std::size_t jp = index(i, wrap(static_cast<long long>(j) + 1, ny_));
            const std::size_t jm = index(i, wrap(static_cast<long long>(j) - 1, ny_));

            const double dudx = (u_[ip] - u_[im]) / (2.0 * dx_);
            const double dudy = (u_[jp] - u_[jm]) / (2.0 * dy_);
            const double d2udx2 = (u_[ip] - 2.0 * u_[idx] + u_[im]) / (dx_ * dx_);
            const double d2udy2 = (u_[jp] - 2.0 * u_[idx] + u_[jm]) / (dy_ * dy_);

            const double dvdx = (v_[ip] - v_[im]) / (2.0 * dx_);
            const double dvdy = (v_[jp] - v_[jm]) / (2.0 * dy_);
            const double d2vdx2 = (v_[ip] - 2.0 * v_[idx] + v_[im]) / (dx_ * dx_);
            const double d2vdy2 = (v_[jp] - 2.0 * v_[idx] + v_[jm]) / (dy_ * dy_);

            const double uConvection = u_[idx] * dudx + v_[idx] * dudy;
            const double vConvection = u_[idx] * dvdx + v_[idx] * dvdy;

            uStar[idx] = u_[idx] + dt * (-uConvection + viscosity_ * (d2udx2 + d2udy2));
            vStar[idx] = v_[idx] + dt * (-vConvection + viscosity_ * (d2vdx2 + d2vdy2));
        }
    }

    projectToDivergenceFree(uStar, vStar, dt);
    time_ += dt;
}

} // namespace aether::solver
