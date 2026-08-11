#include "aether/solver/KEpsilonLidDrivenCavitySolver2D.hpp"

#include "aether/solver/MixingLengthLidDrivenCavitySolver2D.hpp"

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

KEpsilonLidDrivenCavitySolver2D::KEpsilonLidDrivenCavitySolver2D(std::size_t nx, std::size_t ny,
                                                                  double lengthX, double lengthY,
                                                                  double viscosity, double lidVelocity)
    : nx_(nx), ny_(ny), lengthX_(lengthX), lengthY_(lengthY), dx_(lengthX / static_cast<double>(nx)),
      dy_(lengthY / static_cast<double>(ny)), viscosity_(viscosity), lidVelocity_(lidVelocity),
      u_(nx * ny, 0.0), v_(nx * ny, 0.0), p_(nx * ny, 0.0), k_(nx * ny, 0.0), epsilon_(nx * ny, 0.0),
      nut_(nx * ny, 0.0), production_(nx * ny, 0.0) {
    // Warm-start velocity from the simpler mixing-length closure -- same
    // cold-start deadlock reason KEpsilonChannelFlowSolver1D warm-starts
    // from MixingLengthChannelFlowSolver1D (u=0 gives zero production, so
    // k/epsilon would collapse before the velocity field ever develops).
    MixingLengthLidDrivenCavitySolver2D primer(nx, ny, lengthX, lengthY, viscosity, lidVelocity);
    const double primerDt = 0.3 * primer.stableTimeStep();
    for (int s = 0; s < 400; ++s) {
        primer.step(primerDt);
    }
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            u_[index(i, j)] = primer.u(i, j);
            v_[index(i, j)] = primer.v(i, j);
        }
    }

    // Small, uniform initial k/epsilon -- not from external data, just a
    // self-consistent starting guess (epsilon0 from k0 via the same
    // mixing-length-scale relation eps ~ Cmu^0.75 * k^1.5 / L used
    // throughout turbulence modeling to relate a turbulence length scale
    // to k and epsilon).
    const double k0 = 0.01 * lidVelocity * lidVelocity;
    const double lengthScale = 0.1 * std::min(lengthX_, lengthY_);
    const double eps0 = std::pow(kCmu, 0.75) * std::pow(k0, 1.5) / lengthScale;
    std::fill(k_.begin(), k_.end(), k0);
    std::fill(epsilon_.begin(), epsilon_.end(), eps0);
}

double KEpsilonLidDrivenCavitySolver2D::u(std::size_t i, std::size_t j) const { return u_[index(i, j)]; }
double KEpsilonLidDrivenCavitySolver2D::v(std::size_t i, std::size_t j) const { return v_[index(i, j)]; }
double KEpsilonLidDrivenCavitySolver2D::pressure(std::size_t i, std::size_t j) const {
    return p_[index(i, j)];
}
double KEpsilonLidDrivenCavitySolver2D::k(std::size_t i, std::size_t j) const { return k_[index(i, j)]; }
double KEpsilonLidDrivenCavitySolver2D::epsilon(std::size_t i, std::size_t j) const {
    return epsilon_[index(i, j)];
}
double KEpsilonLidDrivenCavitySolver2D::eddyViscosity(std::size_t i, std::size_t j) const {
    return nut_[index(i, j)];
}

double KEpsilonLidDrivenCavitySolver2D::lidVelocityAt(std::size_t i) const {
    constexpr double kPi = 3.14159265358979323846;
    const double x = (static_cast<double>(i) + 0.5) * dx_;
    const double s = std::sin(kPi * x / lengthX_);
    return lidVelocity_ * s * s;
}

double KEpsilonLidDrivenCavitySolver2D::stableTimeStep() const {
    double maxNut = 0.0;
    for (double n : nut_) {
        maxNut = std::max(maxNut, n);
    }
    const double effectiveViscosity = viscosity_ + maxNut;
    const double diffusiveLimit =
        1.0 / (2.0 * effectiveViscosity * (1.0 / (dx_ * dx_) + 1.0 / (dy_ * dy_)));
    const double convectiveLimit = std::min(dx_, dy_) / std::max(lidVelocity_, 1e-12);
    return std::min(diffusiveLimit, convectiveLimit);
}

double KEpsilonLidDrivenCavitySolver2D::dirichletAt(const std::vector<double>& field, std::size_t i,
                                                     std::size_t j, int di, int dj, double wallValue) const {
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

double KEpsilonLidDrivenCavitySolver2D::neumannAt(const std::vector<double>& field, std::size_t i,
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

double KEpsilonLidDrivenCavitySolver2D::epsilonGhostAt(std::size_t i, std::size_t j, int di, int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    const bool crossesWall = ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 ||
                              nj >= static_cast<long long>(ny_);
    if (!crossesWall) {
        return epsilon_[index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
    }
    // The interior cell (i,j) itself sits half a cell from the wall being
    // crossed; low-Re asymptotic relation epsilon_wall = 2*nu*k/y_half^2.
    const double halfSpacing = (di != 0 ? dx_ : dy_) / 2.0;
    const double kNearWall = std::max(k_[index(i, j)], 0.0);
    const double epsWall = 2.0 * viscosity_ * kNearWall / (halfSpacing * halfSpacing);
    return 2.0 * epsWall - epsilon_[index(i, j)];
}

double KEpsilonLidDrivenCavitySolver2D::nutAt(std::size_t i, std::size_t j, int di, int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0 || ni >= static_cast<long long>(nx_) || nj < 0 || nj >= static_cast<long long>(ny_)) {
        return 0.0;
    }
    return nut_[index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
}

void KEpsilonLidDrivenCavitySolver2D::updateEddyViscosityAndProduction() {
    for (std::size_t j = 0; j < ny_; ++j) {
        for (std::size_t i = 0; i < nx_; ++i) {
            const std::size_t idx = index(i, j);
            const double kSafe = std::max(k_[idx], kFloor);
            const double epsSafe = std::max(epsilon_[idx], kFloor);
            nut_[idx] = kCmu * kSafe * kSafe / epsSafe;

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
            const double strainSquared = 2.0 * dudx * dudx + 2.0 * dvdy * dvdy + shear * shear;

            production_[idx] = nut_[idx] * strainSquared;
        }
    }
}

double KEpsilonLidDrivenCavitySolver2D::maxDivergence() const {
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

std::vector<double> KEpsilonLidDrivenCavitySolver2D::applyLaplacian(const std::vector<double>& x) const {
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

double KEpsilonLidDrivenCavitySolver2D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void KEpsilonLidDrivenCavitySolver2D::projectToDivergenceFree(std::vector<double>& uStar,
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
            const double gdpdx = (pE - pW) / (2.0 * dx_);
            const double gdpdy = (pN - pS) / (2.0 * dy_);
            u_[idx] = uStar[idx] - dt * gdpdx;
            v_[idx] = vStar[idx] - dt * gdpdy;
        }
    }
}

void KEpsilonLidDrivenCavitySolver2D::step(double dt) {
    updateEddyViscosityAndProduction();

    std::vector<double> uStar(u_.size());
    std::vector<double> vStar(v_.size());
    std::vector<double> kStar(k_.size());
    std::vector<double> epsStar(epsilon_.size());

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

            // --- k transport ---
            {
                const double kE = dirichletAt(k_, i, j, 1, 0, 0.0);
                const double kW = dirichletAt(k_, i, j, -1, 0, 0.0);
                const double kN = dirichletAt(k_, i, j, 0, 1, 0.0);
                const double kS = dirichletAt(k_, i, j, 0, -1, 0.0);
                const double dkdx = (kE - kW) / (2.0 * dx_);
                const double dkdy = (kN - kS) / (2.0 * dy_);
                const double advectionK = u_[idx] * dkdx + v_[idx] * dkdy;

                const double gEk = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 1, 0)) / kSigmaK;
                const double gWk = viscosity_ + 0.5 * (nutHere + nutAt(i, j, -1, 0)) / kSigmaK;
                const double gNk = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 0, 1)) / kSigmaK;
                const double gSk = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 0, -1)) / kSigmaK;
                const double diffusionK = (gEk * (kE - k_[idx]) - gWk * (k_[idx] - kW)) / (dx_ * dx_) +
                                           (gNk * (kN - k_[idx]) - gSk * (k_[idx] - kS)) / (dy_ * dy_);

                const double newK = k_[idx] + dt * (-advectionK + production_[idx] - epsilon_[idx] +
                                                      diffusionK);
                kStar[idx] = std::max(newK, kFloor);
            }

            // --- epsilon transport (Patankar-style implicit destruction
            // even within this explicit march: the update below solves
            // newEps = numerator/denominator with the quadratic sink
            // -Ceps2*eps^2/k placed in the denominator using the *current*
            // eps/k as the frozen linearization factor, rather than adding
            // an explicit -Ceps2*eps^2/k*dt term directly -- the same
            // stability fix KEpsilonChannelFlowSolver1D needed for its own
            // epsilon destruction term, applied here to avoid re-hitting
            // the identical instability in a transient setting.) ---
            {
                const double epsE = epsilonGhostAt(i, j, 1, 0);
                const double epsW = epsilonGhostAt(i, j, -1, 0);
                const double epsN = epsilonGhostAt(i, j, 0, 1);
                const double epsS = epsilonGhostAt(i, j, 0, -1);
                const double depsdx = (epsE - epsW) / (2.0 * dx_);
                const double depsdy = (epsN - epsS) / (2.0 * dy_);
                const double advectionEps = u_[idx] * depsdx + v_[idx] * depsdy;

                const double gEe = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 1, 0)) / kSigmaEps;
                const double gWe = viscosity_ + 0.5 * (nutHere + nutAt(i, j, -1, 0)) / kSigmaEps;
                const double gNe = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 0, 1)) / kSigmaEps;
                const double gSe = viscosity_ + 0.5 * (nutHere + nutAt(i, j, 0, -1)) / kSigmaEps;
                const double diffusionEps =
                    (gEe * (epsE - epsilon_[idx]) - gWe * (epsilon_[idx] - epsW)) / (dx_ * dx_) +
                    (gNe * (epsN - epsilon_[idx]) - gSe * (epsilon_[idx] - epsS)) / (dy_ * dy_);

                const double kSafe = std::max(k_[idx], kFloor);
                const double productionEps = kCeps1 * (epsilon_[idx] / kSafe) * production_[idx];
                const double destructionRate = kCeps2 * epsilon_[idx] / kSafe; // frozen factor

                const double explicitPart =
                    epsilon_[idx] + dt * (-advectionEps + productionEps + diffusionEps);
                const double newEps = explicitPart / (1.0 + dt * destructionRate);
                epsStar[idx] = std::max(newEps, kFloor);
            }
        }
    }

    k_ = std::move(kStar);
    epsilon_ = std::move(epsStar);

    projectToDivergenceFree(uStar, vStar, dt);
    time_ += dt;
}

} // namespace aether::solver
