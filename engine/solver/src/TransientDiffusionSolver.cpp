#include "aether/solver/TransientDiffusionSolver.hpp"

#include <limits>
#include <utility>

namespace aether::solver {

using aether::core::Vector3;

void TransientDiffusionSolver::setValue(std::size_t i, std::size_t j, std::size_t k, double newValue) {
    field_[grid_->cellIndex(i, j, k)] = newValue;
}

double TransientDiffusionSolver::stableTimeStep() const {
    const Vector3 h = grid_->spacing();
    double totalWeight = 0.0;
    if (grid_->nx() > 1) {
        totalWeight += 1.0 / (h.x * h.x);
    }
    if (grid_->ny() > 1) {
        totalWeight += 1.0 / (h.y * h.y);
    }
    if (grid_->nz() > 1) {
        totalWeight += 1.0 / (h.z * h.z);
    }
    return totalWeight > 0.0 ? 1.0 / totalWeight : std::numeric_limits<double>::infinity();
}

void TransientDiffusionSolver::step(double dt) {
    const auto directions = neighborDirections();

    const std::size_t nx = grid_->nx();
    const std::size_t ny = grid_->ny();
    const std::size_t nz = grid_->nz();

    // A separate buffer is required (not in-place Gauss-Seidel-style
    // updates): explicit time integration must use the *same* time-level
    // values for every neighbor within a step.
    std::vector<double> next = field_;
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

                const double laplacian = weightedSum - weightTotal * field_[idx];
                next[idx] = field_[idx] + dt * (laplacian + source_);
            }
        }
    }
    field_ = std::move(next);
    time_ += dt;
}

} // namespace aether::solver
