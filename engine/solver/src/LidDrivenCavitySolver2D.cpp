#include "aether/solver/LidDrivenCavitySolver2D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

LidDrivenCavitySolver2D::LidDrivenCavitySolver2D(std::size_t nx, std::size_t ny, double lengthX,
                                                   double lengthY, double viscosity, double lidVelocity)
    : nx_(nx), ny_(ny), dx_(lengthX / static_cast<double>(nx)), dy_(lengthY / static_cast<double>(ny)),
      viscosity_(viscosity), lidVelocity_(lidVelocity), u_(nx * ny, 0.0), v_(nx * ny, 0.0),
      p_(nx * ny, 0.0) {}

double LidDrivenCavitySolver2D::u(std::size_t i, std::size_t j) const { return u_[index(i, j)]; }
double LidDrivenCavitySolver2D::v(std::size_t i, std::size_t j) const { return v_[index(i, j)]; }
double LidDrivenCavitySolver2D::pressure(std::size_t i, std::size_t j) const { return p_[index(i, j)]; }

double LidDrivenCavitySolver2D::lidVelocityAt(std::size_t i) const {
    constexpr double kPi = 3.14159265358979323846;
    const double lengthX = dx_ * static_cast<double>(nx_);
    const double x = (static_cast<double>(i) + 0.5) * dx_;
    const double s = std::sin(kPi * x / lengthX);
    return lidVelocity_ * s * s;
}

double LidDrivenCavitySolver2D::stableTimeStep() const {
    const double diffusiveLimit =
        1.0 / (2.0 * viscosity_ * (1.0 / (dx_ * dx_) + 1.0 / (dy_ * dy_)));
    const double convectiveLimit = std::min(dx_, dy_) / std::max(lidVelocity_, 1e-12);
    return std::min(diffusiveLimit, convectiveLimit);
}

double LidDrivenCavitySolver2D::dirichletAt(const std::vector<double>& field, std::size_t i, std::size_t j,
                                              int di, int dj, double wallValue) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0) {
        return -field[index(0, j)]; // left wall: value = 0
    }
    if (ni >= static_cast<long long>(nx_)) {
        return -field[index(nx_ - 1, j)]; // right wall: value = 0
    }
    if (nj < 0) {
        return -field[index(i, 0)]; // bottom wall: value = 0
    }
    if (nj >= static_cast<long long>(ny_)) {
        return 2.0 * wallValue - field[index(i, ny_ - 1)]; // top wall (the lid, for u)
    }
    return field[index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
}

double LidDrivenCavitySolver2D::neumannAt(const std::vector<double>& field, std::size_t i, std::size_t j,
                                            int di, int dj) const {
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

double LidDrivenCavitySolver2D::maxDivergence() const {
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

std::vector<double> LidDrivenCavitySolver2D::applyLaplacian(const std::vector<double>& x) const {
    const double ax = 1.0 / (dx_ * dx_);
    const double ay = 1.0 / (dy_ * dy_);
    const double weightTotal = 2.0 * ax + 2.0 * ay;

    std::vector<double> result(x.size());
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            if (idx == 0) {
                result[idx] = x[idx]; // pinned reference cell: removes the Neumann null space
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

double LidDrivenCavitySolver2D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void LidDrivenCavitySolver2D::projectToDivergenceFree(std::vector<double>& uStar,
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
            // applyLaplacian() computes -nabla^2(p), so the right-hand side
            // of nabla^2(p) = div(uStar)/dt must be negated here (see the
            // sign bug this exact mistake caused in TaylorGreenVortexSolver2D).
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
            const double dpdx = (pE - pW) / (2.0 * dx_);
            const double dpdy = (pN - pS) / (2.0 * dy_);
            u_[idx] = uStar[idx] - dt * dpdx;
            v_[idx] = vStar[idx] - dt * dpdy;
        }
    }
}

void LidDrivenCavitySolver2D::step(double dt) {
    std::vector<double> uStar(u_.size());
    std::vector<double> vStar(v_.size());

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
            const double d2udx2 = (uE - 2.0 * u_[idx] + uW) / (dx_ * dx_);
            const double d2udy2 = (uN - 2.0 * u_[idx] + uS) / (dy_ * dy_);

            const double dvdx = (vE - vW) / (2.0 * dx_);
            const double dvdy = (vN - vS) / (2.0 * dy_);
            const double d2vdx2 = (vE - 2.0 * v_[idx] + vW) / (dx_ * dx_);
            const double d2vdy2 = (vN - 2.0 * v_[idx] + vS) / (dy_ * dy_);

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
