#include "aether/gpu/PoissonOperatorCuda.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using aether::gpu::PoissonOperatorCuda;

namespace {

// Independent CPU reference, deliberately reimplemented here rather than
// calling into the engine (StaggeredCavityBase3D::applyLaplacian is
// protected, not part of its public surface, and re-deriving the formula
// from the class's own documented behavior -- rather than reaching for a
// friend declaration or making it public just for this test -- is the same
// "independent reimplementation as cross-check" discipline this project's
// test suite already uses elsewhere (e.g. satisfiesDelaunayProperty()'s
// own brute-force checks). Any accidental drift between the two would show
// up as a real mismatch here, not get silently hidden by sharing code.
std::vector<double> cpuReferencePoisson(const std::vector<double>& x, std::size_t nx, std::size_t ny,
                                         std::size_t nz, double dx, double dy, double dz) {
    const auto index = [nx, ny](std::size_t i, std::size_t j, std::size_t k) { return i + j * nx + k * nx * ny; };
    const double ax = 1.0 / (dx * dx);
    const double ay = 1.0 / (dy * dy);
    const double az = 1.0 / (dz * dz);
    const double weightTotal = 2.0 * ax + 2.0 * ay + 2.0 * az;

    std::vector<double> result(x.size());
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const std::size_t idx = index(i, j, k);
                if (idx == 0) {
                    result[idx] = x[idx];
                    continue;
                }
                const double left = i > 0 ? x[index(i - 1, j, k)] : x[index(0, j, k)];
                const double right = i + 1 < nx ? x[index(i + 1, j, k)] : x[index(nx - 1, j, k)];
                const double down = j > 0 ? x[index(i, j - 1, k)] : x[index(i, 0, k)];
                const double up = j + 1 < ny ? x[index(i, j + 1, k)] : x[index(i, ny - 1, k)];
                const double back = k > 0 ? x[index(i, j, k - 1)] : x[index(i, j, 0)];
                const double front = k + 1 < nz ? x[index(i, j, k + 1)] : x[index(i, j, nz - 1)];
                const double weightedSum = ax * (left + right) + ay * (down + up) + az * (back + front);
                result[idx] = weightTotal * x[idx] - weightedSum;
            }
        }
    }
    return result;
}

// A deterministic, non-trivial field (no two adjacent cells equal, no
// obvious symmetry that could hide an indexing bug), built only from sin/
// cos/+/-/* so it is reproducible without a random-number generator.
std::vector<double> testField(std::size_t nx, std::size_t ny, std::size_t nz) {
    std::vector<double> field(nx * ny * nz);
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const double x = static_cast<double>(i);
                const double y = static_cast<double>(j);
                const double z = static_cast<double>(k);
                field[i + j * nx + k * nx * ny] =
                    std::sin(0.7 * x + 0.3 * y) + std::cos(0.5 * y - 0.2 * z) + 0.1 * x * z - 0.05 * y * y;
            }
        }
    }
    return field;
}

// **The point of this test**: MSVC and nvcc are both IEEE-754-strict by
// default, and this kernel performs the identical sequence of operations
// the CPU reference does, in the identical order, on identical binary64
// doubles -- so GPU and CPU outputs are expected to agree exactly, not
// just to some tolerance. A nonzero difference here means the two
// implementations are doing genuinely different arithmetic (wrong
// neighbour index, wrong operation order, or accidentally running at
// single precision), not floating-point noise to be tolerated.
void testPoissonOperatorCudaMatchesCpuExactly() {
    const std::size_t nx = 6;
    const std::size_t ny = 5;
    const std::size_t nz = 4;
    const double dx = 0.1;
    const double dy = 0.2;
    const double dz = 0.05; // deliberately anisotropic, exercises ax != ay != az

    PoissonOperatorCuda op(nx, ny, nz, dx, dy, dz);
    if (!op.available()) {
        // The CUDA Toolkit must be present at build time for this
        // executable to exist at all (see engine/gpu/CMakeLists.txt), but
        // that doesn't guarantee a physical device at *run* time -- a
        // headless CI machine could compile this without one attached.
        // Reported explicitly rather than silently passing or hard-failing
        // on hardware this project's own dev machine doesn't need to
        // assume everyone has.
        std::printf("  [aether_gpu_tests] nenhum dispositivo CUDA disponivel em tempo de execucao -- "
                     "pulando as checagens numericas.\n");
        return;
    }

    const std::vector<double> field = testField(nx, ny, nz);
    const std::vector<double> gpuResult = op.apply(field);
    const std::vector<double> cpuResult = cpuReferencePoisson(field, nx, ny, nz, dx, dy, dz);

    AETHER_CHECK(gpuResult.size() == cpuResult.size());
    double maxDiff = 0.0;
    std::size_t mismatches = 0;
    for (std::size_t i = 0; i < gpuResult.size(); ++i) {
        const double diff = std::fabs(gpuResult[i] - cpuResult[i]);
        maxDiff = std::max(maxDiff, diff);
        if (gpuResult[i] != cpuResult[i]) {
            ++mismatches;
        }
    }
    std::printf("  [aether_gpu_tests] GPU vs CPU: %zu/%zu celulas identicas bit a bit, maxDiff=%.3e\n",
                gpuResult.size() - mismatches, gpuResult.size(), maxDiff);

    // Exact equality, not a tolerance -- see the function comment above.
    // (Measured before -fmad=false was added to engine/gpu/CMakeLists.txt:
    // nvcc's default FMA contraction gave ~1e-13..1e-14 absolute
    // differences here, the signature of one extra/missing rounding step,
    // not an algorithmic bug -- see that CMakeLists.txt comment.)
    for (std::size_t i = 0; i < gpuResult.size(); ++i) {
        AETHER_CHECK(gpuResult[i] == cpuResult[i]);
    }

    // Sanity check the reference itself isn't a no-op: the operator must
    // actually mix neighbours, so it should not just return the input
    // back unchanged (a bug that bit-exact equality between two identical
    // no-ops would not catch).
    bool anyDifferentFromInput = false;
    for (std::size_t i = 0; i < field.size(); ++i) {
        if (cpuResult[i] != field[i]) {
            anyDifferentFromInput = true;
            break;
        }
    }
    AETHER_CHECK(anyDifferentFromInput);
}

// The pinned cell (index 0) must be the identity on both sides, for any
// input -- checked separately since testPoissonOperatorCudaMatchesCpuExactly()'s
// field is unlikely to happen to make that visible on its own.
void testPoissonOperatorCudaPinnedCellIsIdentity() {
    const std::size_t nx = 4;
    const std::size_t ny = 4;
    const std::size_t nz = 4;
    PoissonOperatorCuda op(nx, ny, nz, 1.0, 1.0, 1.0);
    if (!op.available()) {
        return; // already reported by the test above
    }
    std::vector<double> field(nx * ny * nz, 3.0);
    field[0] = 7.5;
    const std::vector<double> result = op.apply(field);
    AETHER_CHECK(result[0] == 7.5);
}

} // namespace

int main() {
    testPoissonOperatorCudaMatchesCpuExactly();
    testPoissonOperatorCudaPinnedCellIsIdentity();
    std::printf("aether_gpu_tests: OK\n");
    return 0;
}
