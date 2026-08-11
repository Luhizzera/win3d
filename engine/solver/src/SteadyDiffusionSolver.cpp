#include "aether/solver/SteadyDiffusionSolver.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

std::size_t SteadyDiffusionSolver::solve(std::size_t maxIterations, double tolerance) {
    const auto directions = neighborDirections();

    const std::size_t nx = grid_->nx();
    const std::size_t ny = grid_->ny();
    const std::size_t nz = grid_->nz();

    std::size_t iteration = 0;
    for (; iteration < maxIterations; ++iteration) {
        double maxDelta = 0.0;
        for (std::size_t k = 0; k < nz; ++k) {
            for (std::size_t j = 0; j < ny; ++j) {
                for (std::size_t i = 0; i < nx; ++i) {
                    const std::size_t idx = grid_->cellIndex(i, j, k);
                    if (isFixed_[idx]) {
                        continue;
                    }

                    double weightedSum = 0.0;
                    double weightTotal = 0.0;
                    for (const auto& [d, weight] : directions) {
                        if (grid_->hasNeighbor(i, j, k, d[0], d[1], d[2])) {
                            const auto ni = static_cast<std::size_t>(static_cast<long long>(i) + d[0]);
                            const auto nj = static_cast<std::size_t>(static_cast<long long>(j) + d[1]);
                            const auto nk = static_cast<std::size_t>(static_cast<long long>(k) + d[2]);
                            weightedSum += weight * field_[grid_->cellIndex(ni, nj, nk)];
                            weightTotal += weight;
                        }
                    }
                    if (weightTotal == 0.0) {
                        continue;
                    }

                    const double newValue = (weightedSum + source_) / weightTotal;
                    maxDelta = std::max(maxDelta, std::fabs(newValue - field_[idx]));
                    field_[idx] = newValue;
                }
            }
        }
        if (maxDelta < tolerance) {
            ++iteration;
            break;
        }
    }
    return iteration;
}

std::vector<double> SteadyDiffusionSolver::applyOperator(const std::vector<double>& x) const {
    const auto directions = neighborDirections();

    const std::size_t nx = grid_->nx();
    const std::size_t ny = grid_->ny();
    const std::size_t nz = grid_->nz();

    std::vector<double> result(x.size(), 0.0);
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t idx = grid_->cellIndex(i, j, k);
                if (isFixed_[idx]) {
                    result[idx] = x[idx]; // identity row: fixed cells assert x == prescribed value
                    continue;
                }

                double weightedSum = 0.0;
                double weightTotal = 0.0;
                for (const auto& [d, weight] : directions) {
                    if (grid_->hasNeighbor(i, j, k, d[0], d[1], d[2])) {
                        const auto ni = static_cast<std::size_t>(static_cast<long long>(i) + d[0]);
                        const auto nj = static_cast<std::size_t>(static_cast<long long>(j) + d[1]);
                        const auto nk = static_cast<std::size_t>(static_cast<long long>(k) + d[2]);
                        weightedSum += weight * x[grid_->cellIndex(ni, nj, nk)];
                        weightTotal += weight;
                    }
                }
                result[idx] = weightTotal * x[idx] - weightedSum;
            }
        }
    }
    return result;
}

std::vector<double> SteadyDiffusionSolver::computeDiagonal() const {
    const auto directions = neighborDirections();

    const std::size_t nx = grid_->nx();
    const std::size_t ny = grid_->ny();
    const std::size_t nz = grid_->nz();

    std::vector<double> diagonal(field_.size(), 1.0);
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t idx = grid_->cellIndex(i, j, k);
                if (isFixed_[idx]) {
                    continue; // identity row: diagonal already 1.0
                }
                double weightTotal = 0.0;
                for (const auto& [d, weight] : directions) {
                    if (grid_->hasNeighbor(i, j, k, d[0], d[1], d[2])) {
                        weightTotal += weight;
                    }
                }
                diagonal[idx] = weightTotal;
            }
        }
    }
    return diagonal;
}

double SteadyDiffusionSolver::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

std::size_t SteadyDiffusionSolver::solveConjugateGradient(std::size_t maxIterations, double tolerance) {
    const std::size_t n = field_.size();

    std::vector<double> rhs(n);
    for (std::size_t idx = 0; idx < n; ++idx) {
        rhs[idx] = isFixed_[idx] ? field_[idx] : source_;
    }

    std::vector<double> residual(n);
    {
        const std::vector<double> operatorAppliedToField = applyOperator(field_);
        for (std::size_t idx = 0; idx < n; ++idx) {
            residual[idx] = isFixed_[idx] ? 0.0 : rhs[idx] - operatorAppliedToField[idx];
        }
    }
    std::vector<double> direction = residual;
    double residualDotResidual = dot(residual, residual);

    std::size_t iteration = 0;
    for (; iteration < maxIterations; ++iteration) {
        if (std::sqrt(residualDotResidual) < tolerance) {
            break;
        }

        const std::vector<double> operatorAppliedToDirection = applyOperator(direction);
        const double directionDotOperatorDirection = dot(direction, operatorAppliedToDirection);
        if (directionDotOperatorDirection == 0.0) {
            break; // degenerate system (e.g. no free cells at all)
        }
        const double alpha = residualDotResidual / directionDotOperatorDirection;

        for (std::size_t idx = 0; idx < n; ++idx) {
            if (!isFixed_[idx]) {
                field_[idx] += alpha * direction[idx];
            }
        }

        std::vector<double> newResidual(n);
        for (std::size_t idx = 0; idx < n; ++idx) {
            newResidual[idx] = isFixed_[idx] ? 0.0 : residual[idx] - alpha * operatorAppliedToDirection[idx];
        }
        const double newResidualDotResidual = dot(newResidual, newResidual);
        const double beta = newResidualDotResidual / residualDotResidual;

        for (std::size_t idx = 0; idx < n; ++idx) {
            direction[idx] = isFixed_[idx] ? 0.0 : newResidual[idx] + beta * direction[idx];
        }

        residual = std::move(newResidual);
        residualDotResidual = newResidualDotResidual;
    }
    return iteration;
}

std::size_t SteadyDiffusionSolver::solvePreconditionedConjugateGradient(std::size_t maxIterations,
                                                                          double tolerance) {
    const std::size_t n = field_.size();
    const std::vector<double> diagonal = computeDiagonal();

    std::vector<double> rhs(n);
    for (std::size_t idx = 0; idx < n; ++idx) {
        rhs[idx] = isFixed_[idx] ? field_[idx] : source_;
    }

    std::vector<double> residual(n);
    {
        const std::vector<double> operatorAppliedToField = applyOperator(field_);
        for (std::size_t idx = 0; idx < n; ++idx) {
            residual[idx] = isFixed_[idx] ? 0.0 : rhs[idx] - operatorAppliedToField[idx];
        }
    }

    std::vector<double> z(n); // preconditioned residual: M^-1 * residual, M = diag(A)
    for (std::size_t idx = 0; idx < n; ++idx) {
        z[idx] = isFixed_[idx] ? 0.0 : residual[idx] / diagonal[idx];
    }
    std::vector<double> direction = z;
    double residualDotZ = dot(residual, z);

    std::size_t iteration = 0;
    for (; iteration < maxIterations; ++iteration) {
        if (std::sqrt(dot(residual, residual)) < tolerance) {
            break;
        }

        const std::vector<double> operatorAppliedToDirection = applyOperator(direction);
        const double directionDotOperatorDirection = dot(direction, operatorAppliedToDirection);
        if (directionDotOperatorDirection == 0.0) {
            break; // degenerate system (e.g. no free cells at all)
        }
        const double alpha = residualDotZ / directionDotOperatorDirection;

        for (std::size_t idx = 0; idx < n; ++idx) {
            if (!isFixed_[idx]) {
                field_[idx] += alpha * direction[idx];
            }
        }

        std::vector<double> newResidual(n);
        for (std::size_t idx = 0; idx < n; ++idx) {
            newResidual[idx] = isFixed_[idx] ? 0.0 : residual[idx] - alpha * operatorAppliedToDirection[idx];
        }

        std::vector<double> newZ(n);
        for (std::size_t idx = 0; idx < n; ++idx) {
            newZ[idx] = isFixed_[idx] ? 0.0 : newResidual[idx] / diagonal[idx];
        }
        const double newResidualDotZ = dot(newResidual, newZ);
        const double beta = newResidualDotZ / residualDotZ;

        for (std::size_t idx = 0; idx < n; ++idx) {
            direction[idx] = isFixed_[idx] ? 0.0 : newZ[idx] + beta * direction[idx];
        }

        residual = std::move(newResidual);
        z = std::move(newZ);
        residualDotZ = newResidualDotZ;
    }
    return iteration;
}

} // namespace aether::solver
