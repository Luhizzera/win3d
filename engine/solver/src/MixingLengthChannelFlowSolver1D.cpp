#include "aether/solver/MixingLengthChannelFlowSolver1D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

namespace {
constexpr double kKarman = 0.41;          // von Karman constant
constexpr double kLengthCapFactor = 0.09; // Escudier's asymptotic mixing-length cap, x half-height
} // namespace

MixingLengthChannelFlowSolver1D::MixingLengthChannelFlowSolver1D(std::size_t ny, double height,
                                                                    double kinematicViscosity,
                                                                    double source)
    : ny_(ny), height_(height), h_(height / static_cast<double>(ny)), nu_(kinematicViscosity),
      source_(source), u_(ny, 0.0), nut_(ny, 0.0) {}

double MixingLengthChannelFlowSolver1D::cellCenterY(std::size_t j) const {
    return (static_cast<double>(j) + 0.5) * h_;
}

double MixingLengthChannelFlowSolver1D::wallDistance(std::size_t j) const {
    const double y = cellCenterY(j);
    return std::min(y, height_ - y);
}

double MixingLengthChannelFlowSolver1D::uAt(long long j) const {
    if (j < 0) {
        return -u_[0]; // Dirichlet mirror: u = 0 at the wall below cell 0
    }
    if (j >= static_cast<long long>(ny_)) {
        return -u_[ny_ - 1]; // Dirichlet mirror: u = 0 at the wall above cell ny_-1
    }
    return u_[static_cast<std::size_t>(j)];
}

double MixingLengthChannelFlowSolver1D::velocityGradientAt(std::size_t j) const {
    return (uAt(static_cast<long long>(j) + 1) - uAt(static_cast<long long>(j) - 1)) / (2.0 * h_);
}

double MixingLengthChannelFlowSolver1D::gammaAt(long long j) const {
    // The wall itself has zero eddy viscosity (mixing length -> 0 there
    // by definition), regardless of the adjacent cell's own nu_t.
    if (j < 0 || j >= static_cast<long long>(ny_)) {
        return nu_;
    }
    return nu_ + nut_[static_cast<std::size_t>(j)];
}

void MixingLengthChannelFlowSolver1D::updateEddyViscosity() {
    const double halfHeight = height_ / 2.0;
    for (std::size_t j = 0; j < ny_; ++j) {
        const double mixingLength = std::min(kKarman * wallDistance(j), kLengthCapFactor * halfHeight);
        const double gradient = velocityGradientAt(j);
        nut_[j] = mixingLength * mixingLength * std::fabs(gradient);
    }
}

std::size_t MixingLengthChannelFlowSolver1D::solve(std::size_t maxOuterIterations, double tolerance) {
    std::size_t outer = 0;
    for (; outer < maxOuterIterations; ++outer) {
        updateEddyViscosity();

        double maxDelta = 0.0;
        // A handful of Gauss-Seidel sweeps per Picard (outer) iteration is
        // enough since nu_t only changes slowly between outer iterations;
        // re-deriving it every single inner sweep would also converge, just
        // less predictably.
        for (int inner = 0; inner < 50; ++inner) {
            maxDelta = 0.0;
            for (std::size_t j = 0; j < ny_; ++j) {
                const double gammaHere = nu_ + nut_[j];
                const double gammaUp = 0.5 * (gammaHere + gammaAt(static_cast<long long>(j) + 1));
                const double gammaDown = 0.5 * (gammaHere + gammaAt(static_cast<long long>(j) - 1));
                const double uUp = uAt(static_cast<long long>(j) + 1);
                const double uDown = uAt(static_cast<long long>(j) - 1);
                const double newU =
                    (gammaUp * uUp + gammaDown * uDown + source_ * h_ * h_) / (gammaUp + gammaDown);
                maxDelta = std::max(maxDelta, std::fabs(newU - u_[j]));
                u_[j] = newU;
            }
        }

        if (maxDelta < tolerance) {
            ++outer;
            break;
        }
    }
    return outer;
}

double MixingLengthChannelFlowSolver1D::frictionVelocity() const {
    return std::sqrt(source_ * height_ / 2.0);
}

} // namespace aether::solver
