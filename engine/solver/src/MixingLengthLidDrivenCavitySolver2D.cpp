#include "aether/solver/ExplicitTimeStep.hpp"
#include "aether/solver/MixingLengthLidDrivenCavitySolver2D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

namespace {
constexpr double kKarman = 0.41;
constexpr double kEscudierFactor = 0.09;
} // namespace

MixingLengthLidDrivenCavitySolver2D::MixingLengthLidDrivenCavitySolver2D(std::size_t nx, std::size_t ny,
                                                                          double lengthX, double lengthY,
                                                                          double viscosity,
                                                                          double lidVelocity)
    : nx_(nx), ny_(ny), lengthX_(lengthX), lengthY_(lengthY), dx_(lengthX / static_cast<double>(nx)),
      dy_(lengthY / static_cast<double>(ny)), viscosity_(viscosity), lidVelocity_(lidVelocity),
      u_(nx * ny, 0.0), v_(nx * ny, 0.0), p_(nx * ny, 0.0), nut_(nx * ny, 0.0) {}

double MixingLengthLidDrivenCavitySolver2D::u(std::size_t i, std::size_t j) const { return u_[index(i, j)]; }
double MixingLengthLidDrivenCavitySolver2D::v(std::size_t i, std::size_t j) const { return v_[index(i, j)]; }
double MixingLengthLidDrivenCavitySolver2D::pressure(std::size_t i, std::size_t j) const {
    return p_[index(i, j)];
}
double MixingLengthLidDrivenCavitySolver2D::eddyViscosity(std::size_t i, std::size_t j) const {
    return nut_[index(i, j)];
}

double MixingLengthLidDrivenCavitySolver2D::lidVelocityAt(std::size_t i) const {
    constexpr double kPi = 3.14159265358979323846;
    const double x = (static_cast<double>(i) + 0.5) * dx_;
    const double s = std::sin(kPi * x / lengthX_);
    return lidVelocity_ * s * s;
}

double MixingLengthLidDrivenCavitySolver2D::wallDistanceAt(std::size_t i, std::size_t j) const {
    const double x = (static_cast<double>(i) + 0.5) * dx_;
    const double y = (static_cast<double>(j) + 0.5) * dy_;
    return std::min({x, lengthX_ - x, y, lengthY_ - y});
}

double MixingLengthLidDrivenCavitySolver2D::stableTimeStep() const {
    double maxNut = 0.0;
    for (double n : nut_) {
        maxNut = std::max(maxNut, n);
    }
    const double effectiveViscosity = viscosity_ + maxNut;
    return explicitStableTimeStep(effectiveViscosity, lidVelocity_, {dx_, dy_});
}

double MixingLengthLidDrivenCavitySolver2D::dirichletAt(const std::vector<double>& field, std::size_t i,
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

double MixingLengthLidDrivenCavitySolver2D::neumannAt(const std::vector<double>& field, std::size_t i,
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

double MixingLengthLidDrivenCavitySolver2D::nutAt(std::size_t i, std::size_t j, int di, int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 || nj >= static_cast<long long>(ny_)) {
        return 0.0; // nu_t vanishes exactly at a solid wall
    }
    return nut_[index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
}

void MixingLengthLidDrivenCavitySolver2D::updateEddyViscosity() {
    const double lengthCap = kEscudierFactor * std::min(lengthX_, lengthY_) / 2.0;
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
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
            const double strainMagnitude =
                std::sqrt(2.0 * dudx * dudx + 2.0 * dvdy * dvdy + shear * shear);

            const double wallDist = wallDistanceAt(i, j);
            const double mixingLength = std::min(kKarman * wallDist, lengthCap);

            nut_[index(i, j)] = mixingLength * mixingLength * strainMagnitude;
        }
    }
}

double MixingLengthLidDrivenCavitySolver2D::maxDivergence() const {
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

std::vector<double> MixingLengthLidDrivenCavitySolver2D::applyLaplacian(const std::vector<double>& x) const {
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

double MixingLengthLidDrivenCavitySolver2D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void MixingLengthLidDrivenCavitySolver2D::projectToDivergenceFree(std::vector<double>& uStar,
                                                                    std::vector<double>& vStar,
                                                                    double dt) {
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
            const double dpdx = (pE - pW) / (2.0 * dx_);
            const double dpdy = (pN - pS) / (2.0 * dy_);
            u_[idx] = uStar[idx] - dt * dpdx;
            v_[idx] = vStar[idx] - dt * dpdy;
        }
    }
}

void MixingLengthLidDrivenCavitySolver2D::step(double dt) {
    updateEddyViscosity();

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
            const double dvdx = (vE - vW) / (2.0 * dx_);
            const double dvdy = (vN - vS) / (2.0 * dy_);

            const double uConvection = u_[idx] * dudx + v_[idx] * dudy;
            const double vConvection = u_[idx] * dvdx + v_[idx] * dvdy;

            // Variable-coefficient diffusion: (nu + nu_t) face-averaged
            // between this cell and each neighbor (nu_t exactly 0.0 across
            // any wall face, per nutAt()), applied independently per
            // component (see the class comment for the simplification this
            // implies relative to a full variable-viscosity stress tensor).
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
        }
    }

    projectToDivergenceFree(uStar, vStar, dt);
    time_ += dt;
}

} // namespace aether::solver
