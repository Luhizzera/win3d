#include "aether/analysis/FlowDiagnostics.hpp"
#include "aether/solver/LidDrivenCavitySolver2D.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using aether::analysis::checkerboardIndex;
using aether::analysis::computeStatistics;
using aether::analysis::maxCourantNumber;
using aether::analysis::summarizeField;
using aether::solver::LidDrivenCavitySolver2D;

namespace {

bool nearlyEqual(double a, double b, double tol = 1e-9) { return std::fabs(a - b) < tol; }

void testComputeStatisticsExact() {
    const std::vector<double> field = {3.0, -1.0, 4.0, 1.0, 5.0};
    const auto stats = computeStatistics(field);
    AETHER_CHECK(stats.minValue == -1.0);
    AETHER_CHECK(stats.maxValue == 5.0);
    AETHER_CHECK(nearlyEqual(stats.mean, 12.0 / 5.0));
}

void testComputeStatisticsThrowsOnEmpty() {
    bool threw = false;
    try {
        computeStatistics({});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    AETHER_CHECK(threw);
}

void testMaxCourantNumberExact() {
    // Cell 0: dt*(|2|/0.1 + |0|/0.2) = 0.1*20 = 2.0
    // Cell 1: dt*(|0|/0.1 + |1|/0.2) = 0.1*5 = 0.5
    const std::vector<double> u = {2.0, 0.0};
    const std::vector<double> v = {0.0, 1.0};
    const double cfl = maxCourantNumber(u, v, 0.1, 0.2, 0.1);
    AETHER_CHECK(nearlyEqual(cfl, 2.0));
}

// Any locally linear field -- constant or a plain ramp -- has zero
// checkerboard residual at every interior cell by construction (a cell's
// value exactly equals the average of its neighbours when the field is
// linear), so the index must come out to *exactly* 0.0, not just small.
void testCheckerboardIndexIsZeroForSmoothFields() {
    const std::size_t nx = 5;
    const std::size_t ny = 5;

    std::vector<double> constant(nx * ny, 7.5);
    AETHER_CHECK(checkerboardIndex(constant, nx, ny) == 0.0);

    std::vector<double> ramp(nx * ny);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            ramp[i + j * nx] = static_cast<double>(i) + 2.0 * static_cast<double>(j);
        }
    }
    AETHER_CHECK(nearlyEqual(checkerboardIndex(ramp, nx, ny), 0.0, 1e-12));
}

// A perfect single-cell checkerboard field(i,j) = A*(-1)^(i+j): every
// neighbour of an interior cell has the opposite sign, so the residual at
// that cell is exactly 2*field(i,j) -- see the header comment's proof.
// The normalized index should therefore land at 1.0.
void testCheckerboardIndexIsOneForPerfectCheckerboard() {
    const std::size_t nx = 6;
    const std::size_t ny = 6;
    std::vector<double> field(nx * ny);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            field[i + j * nx] = 3.0 * sign;
        }
    }
    const double index = checkerboardIndex(field, nx, ny);
    std::printf("  [aether_analysis_tests] checkerboardIndex(perfect checkerboard) = %.12f\n", index);
    AETHER_CHECK(nearlyEqual(index, 1.0, 1e-9));
}

void testCheckerboardIndexThrowsOnTooSmallGrid() {
    bool threw = false;
    try {
        checkerboardIndex({1.0, 2.0}, 2, 1);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    AETHER_CHECK(threw);
}

void testSummarizeFieldMatchesStatistics() {
    const std::vector<double> field = {1.0, 2.0, 3.0};
    const std::string text = summarizeField("u", field);
    AETHER_CHECK(text.find("u:") != std::string::npos);
    AETHER_CHECK(text.find("min=1") != std::string::npos);
    AETHER_CHECK(text.find("max=3") != std::string::npos);
    AETHER_CHECK(text.find("n=3") != std::string::npos);
}

// The point of this module: measure the collocated-grid checkerboard
// weakness LidDrivenCavitySolver2D's own class comment documents, on a
// real solved pressure field, rather than only on synthetic extremes.
// A simple one-pass 4-neighbour smoothing of the interior cells should
// strictly *reduce* the measured index relative to the raw field --
// smoothing removes exactly the kind of cell-to-cell sign alternation the
// index is designed to detect, so this is a real, checkable claim, not
// just "the number changed".
void testCheckerboardIndexDropsAfterSmoothingRealPressureField() {
    const std::size_t nx = 20;
    const std::size_t ny = 16;
    LidDrivenCavitySolver2D solver(nx, ny, 1.0, 0.8, 0.01, 1.0);
    for (int step = 0; step < 300; ++step) {
        solver.step(solver.stableTimeStep());
    }

    std::vector<double> pressure(nx * ny);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            pressure[i + j * nx] = solver.pressure(i, j);
        }
    }

    std::vector<double> smoothed = pressure;
    for (std::size_t j = 1; j + 1 < ny; ++j) {
        for (std::size_t i = 1; i + 1 < nx; ++i) {
            smoothed[i + j * nx] = 0.25 * (pressure[(i - 1) + j * nx] + pressure[(i + 1) + j * nx] +
                                            pressure[i + (j - 1) * nx] + pressure[i + (j + 1) * nx]);
        }
    }

    const double rawIndex = checkerboardIndex(pressure, nx, ny);
    const double smoothedIndex = checkerboardIndex(smoothed, nx, ny);
    std::printf("  [aether_analysis_tests] real cavity pressure: raw checkerboardIndex=%.6f, "
                "smoothed=%.6f\n",
                rawIndex, smoothedIndex);

    AETHER_CHECK(smoothedIndex < rawIndex);
}

} // namespace

int main() {
    testComputeStatisticsExact();
    testComputeStatisticsThrowsOnEmpty();
    testMaxCourantNumberExact();
    testCheckerboardIndexIsZeroForSmoothFields();
    testCheckerboardIndexIsOneForPerfectCheckerboard();
    testCheckerboardIndexThrowsOnTooSmallGrid();
    testSummarizeFieldMatchesStatistics();
    testCheckerboardIndexDropsAfterSmoothingRealPressureField();
    std::printf("aether_analysis_tests: OK\n");
    return 0;
}
