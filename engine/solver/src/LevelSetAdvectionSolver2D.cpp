#include "aether/solver/LevelSetAdvectionSolver2D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

LevelSetAdvectionSolver2D::LevelSetAdvectionSolver2D(std::size_t nx, std::size_t ny, double lengthX,
                                                       double lengthY)
    : nx_(nx), ny_(ny), lengthX_(lengthX), lengthY_(lengthY), dx_(lengthX / static_cast<double>(nx)),
      dy_(lengthY / static_cast<double>(ny)), phi_(nx * ny, 0.0), u_(nx * ny, 0.0), v_(nx * ny, 0.0) {}

double LevelSetAdvectionSolver2D::cellCenterX(std::size_t i) const {
    return (static_cast<double>(i) + 0.5) * dx_;
}

double LevelSetAdvectionSolver2D::cellCenterY(std::size_t j) const {
    return (static_cast<double>(j) + 0.5) * dy_;
}

void LevelSetAdvectionSolver2D::initialize(const std::function<double(double, double)>& signedDistance) {
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            phi_[index(i, j)] = signedDistance(cellCenterX(i), cellCenterY(j));
        }
    }
}

void LevelSetAdvectionSolver2D::setVelocityField(
    const std::function<void(double, double, double&, double&)>& velocity) {
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            velocity(cellCenterX(i), cellCenterY(j), u_[index(i, j)], v_[index(i, j)]);
        }
    }
}

double LevelSetAdvectionSolver2D::phiAt(long long i, long long j) const {
    const long long clampedI = std::clamp<long long>(i, 0, static_cast<long long>(nx_) - 1);
    const long long clampedJ = std::clamp<long long>(j, 0, static_cast<long long>(ny_) - 1);
    return phi_[index(static_cast<std::size_t>(clampedI), static_cast<std::size_t>(clampedJ))];
}

void LevelSetAdvectionSolver2D::step(double dt) {
    std::vector<double> next(phi_.size());
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const long long ii = static_cast<long long>(i);
            const long long jj = static_cast<long long>(j);
            const double u = u_[index(i, j)];
            const double v = v_[index(i, j)];
            const double phiHere = phiAt(ii, jj);
            const double dphidx =
                (u >= 0.0) ? (phiHere - phiAt(ii - 1, jj)) / dx_ : (phiAt(ii + 1, jj) - phiHere) / dx_;
            const double dphidy =
                (v >= 0.0) ? (phiHere - phiAt(ii, jj - 1)) / dy_ : (phiAt(ii, jj + 1) - phiHere) / dy_;
            next[index(i, j)] = phiHere - dt * (u * dphidx + v * dphidy);
        }
    }
    phi_ = std::move(next);
    time_ += dt;
}

double LevelSetAdvectionSolver2D::stableTimeStep(double cfl) const {
    double maxSpeed = 0.0;
    for (std::size_t k = 0; k < u_.size(); ++k) {
        maxSpeed = std::max({maxSpeed, std::fabs(u_[k]), std::fabs(v_[k])});
    }
    return cfl * std::min(dx_, dy_) / maxSpeed;
}

double LevelSetAdvectionSolver2D::insideArea() const {
    std::size_t insideCells = 0;
    for (double value : phi_) {
        if (value > 0.0) {
            ++insideCells;
        }
    }
    return static_cast<double>(insideCells) * dx_ * dy_;
}

} // namespace aether::solver
