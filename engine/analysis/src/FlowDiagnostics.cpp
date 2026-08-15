#include "aether/analysis/FlowDiagnostics.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace aether::analysis {

FieldStatistics computeStatistics(const std::vector<double>& field) {
    if (field.empty()) {
        throw std::invalid_argument("computeStatistics: field must not be empty");
    }
    const auto [minIt, maxIt] = std::minmax_element(field.begin(), field.end());
    const double sum = std::accumulate(field.begin(), field.end(), 0.0);
    return FieldStatistics{*minIt, *maxIt, sum / static_cast<double>(field.size())};
}

double maxCourantNumber(const std::vector<double>& u, const std::vector<double>& v, double dx, double dy,
                         double dt) {
    if (u.size() != v.size()) {
        throw std::invalid_argument("maxCourantNumber: u and v must be the same size");
    }
    double maxCfl = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        const double cfl = dt * (std::fabs(u[i]) / dx + std::fabs(v[i]) / dy);
        maxCfl = std::max(maxCfl, cfl);
    }
    return maxCfl;
}

double checkerboardIndex(const std::vector<double>& field, std::size_t nx, std::size_t ny) {
    if (field.size() != nx * ny) {
        throw std::invalid_argument("checkerboardIndex: field size does not match nx*ny");
    }
    if (nx < 3 || ny < 3) {
        throw std::invalid_argument("checkerboardIndex: nx and ny must both be >= 3 (need interior cells)");
    }

    const auto index = [nx](std::size_t i, std::size_t j) { return i + j * nx; };

    const double mean = std::accumulate(field.begin(), field.end(), 0.0) / static_cast<double>(field.size());
    double varianceSum = 0.0;
    for (double value : field) {
        const double d = value - mean;
        varianceSum += d * d;
    }
    const double fieldRms = std::sqrt(varianceSum / static_cast<double>(field.size()));
    if (fieldRms == 0.0) {
        return 0.0;
    }

    double residualSquaredSum = 0.0;
    std::size_t interiorCount = 0;
    for (std::size_t j = 1; j + 1 < ny; ++j) {
        for (std::size_t i = 1; i + 1 < nx; ++i) {
            const double neighborAverage =
                0.25 * (field[index(i - 1, j)] + field[index(i + 1, j)] + field[index(i, j - 1)] +
                        field[index(i, j + 1)]);
            const double residual = field[index(i, j)] - neighborAverage;
            residualSquaredSum += residual * residual;
            ++interiorCount;
        }
    }
    const double residualRms = std::sqrt(residualSquaredSum / static_cast<double>(interiorCount));

    return residualRms / (2.0 * fieldRms);
}

std::string summarizeField(const std::string& name, const std::vector<double>& field) {
    const FieldStatistics stats = computeStatistics(field);
    std::ostringstream out;
    out << name << ": min=" << stats.minValue << " max=" << stats.maxValue << " mean=" << stats.mean
        << " (n=" << field.size() << ")";
    return out.str();
}

} // namespace aether::analysis
