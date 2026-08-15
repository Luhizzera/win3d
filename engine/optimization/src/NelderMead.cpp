#include "aether/optimization/NelderMead.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace aether::optimization {

namespace {

std::vector<double> addScaled(const std::vector<double>& a, const std::vector<double>& b, double scaleB) {
    std::vector<double> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + scaleB * b[i];
    }
    return result;
}

std::vector<double> subtract(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

} // namespace

NelderMead::NelderMead(std::size_t maxIterations, double tolerance)
    : maxIterations_(maxIterations), tolerance_(tolerance) {}

OptimizationResult NelderMead::minimize(const ObjectiveFunction& objective, const std::vector<double>& initialGuess,
                                         double initialStepSize) const {
    const std::size_t n = initialGuess.size();
    if (n == 0) {
        throw std::invalid_argument("NelderMead::minimize: initialGuess must be non-empty");
    }

    // Initial simplex: initialGuess plus one vertex perturbed along each
    // axis by initialStepSize.
    std::vector<std::vector<double>> vertices(n + 1, initialGuess);
    for (std::size_t i = 0; i < n; ++i) {
        vertices[i + 1][i] += initialStepSize;
    }
    std::vector<double> values(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        values[i] = objective(vertices[i]);
    }

    constexpr double kAlpha = 1.0; // reflection
    constexpr double kGamma = 2.0; // expansion
    constexpr double kRho = 0.5;   // contraction
    constexpr double kSigma = 0.5; // shrink

    std::vector<std::size_t> order(n + 1);
    const auto sortByValue = [&order, &values]() {
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&values](std::size_t a, std::size_t b) { return values[a] < values[b]; });
    };

    std::size_t iteration = 0;
    bool converged = false;

    for (; iteration < maxIterations_; ++iteration) {
        sortByValue();

        const double best = values[order[0]];
        const double worst = values[order[n]];
        if (std::fabs(worst - best) < tolerance_) {
            converged = true;
            break;
        }

        std::vector<double> centroid(n, 0.0);
        for (std::size_t k = 0; k < n; ++k) {
            const std::vector<double>& v = vertices[order[k]];
            for (std::size_t d = 0; d < n; ++d) {
                centroid[d] += v[d];
            }
        }
        for (double& c : centroid) {
            c /= static_cast<double>(n);
        }

        const std::size_t worstIndex = order[n];
        const std::vector<double> worstVertex = vertices[worstIndex];

        const std::vector<double> reflected = addScaled(centroid, subtract(centroid, worstVertex), kAlpha);
        const double reflectedValue = objective(reflected);

        if (reflectedValue < values[order[0]]) {
            const std::vector<double> expanded = addScaled(centroid, subtract(reflected, centroid), kGamma);
            const double expandedValue = objective(expanded);
            if (expandedValue < reflectedValue) {
                vertices[worstIndex] = expanded;
                values[worstIndex] = expandedValue;
            } else {
                vertices[worstIndex] = reflected;
                values[worstIndex] = reflectedValue;
            }
            continue;
        }

        if (reflectedValue < values[order[n - 1]]) {
            vertices[worstIndex] = reflected;
            values[worstIndex] = reflectedValue;
            continue;
        }

        bool shrink = false;
        if (reflectedValue < worst) {
            // Outside contraction: the reflected point did improve on the
            // worst vertex, just not enough to accept outright.
            const std::vector<double> contracted = addScaled(centroid, subtract(reflected, centroid), kRho);
            const double contractedValue = objective(contracted);
            if (contractedValue <= reflectedValue) {
                vertices[worstIndex] = contracted;
                values[worstIndex] = contractedValue;
            } else {
                shrink = true;
            }
        } else {
            // Inside contraction: even reflecting made things worse.
            const std::vector<double> contracted = addScaled(centroid, subtract(worstVertex, centroid), kRho);
            const double contractedValue = objective(contracted);
            if (contractedValue < worst) {
                vertices[worstIndex] = contracted;
                values[worstIndex] = contractedValue;
            } else {
                shrink = true;
            }
        }

        if (shrink) {
            const std::vector<double> bestVertex = vertices[order[0]];
            for (std::size_t k = 1; k <= n; ++k) {
                const std::size_t idx = order[k];
                vertices[idx] = addScaled(bestVertex, subtract(vertices[idx], bestVertex), kSigma);
                values[idx] = objective(vertices[idx]);
            }
        }
    }

    sortByValue();
    const std::size_t bestIndex = order[0];
    return OptimizationResult{vertices[bestIndex], values[bestIndex], iteration, converged};
}

} // namespace aether::optimization
