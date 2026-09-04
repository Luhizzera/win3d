#include "aether/gpu/ConjugateGradientSolverCuda.hpp"
#include "aether/gpu/MomentumPredictorCuda.hpp"
#include "aether/gpu/PoissonOperatorCuda.hpp"
#include "aether/testing/Check.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using aether::gpu::ConjugateGradientSolverCuda;
using aether::gpu::MomentumPredictorCuda;
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

double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

struct CpuCgResult {
    std::vector<double> pressure;
    std::size_t iterations = 0;
    double residualNorm = 0.0;
    bool brokeDown = false;
};

// Independent CPU reference for the whole CG loop, not just the stencil --
// built on cpuReferencePoisson() above rather than calling into
// aether_solver (same reasoning as that function's own comment: this
// target links only aether_gpu + aether_testing). Structurally the same
// loop as StaggeredCavityBase3D::projectToDivergenceFree
// (engine/solver/src/StaggeredCavityBase3D.cpp), reduced to updating
// `residual`/`direction` in place rather than via a separate
// "newResidual" array -- the two are equivalent (no cell's post-update
// value is needed pre-update anywhere else in the same step), and this is
// the same in-place structure ConjugateGradientSolverCuda's own kernels
// use, which is exactly what this reference exists to cross-check.
CpuCgResult cpuReferenceConjugateGradient(const std::vector<double>& rhs, const std::vector<double>& initialGuess,
                                           std::size_t nx, std::size_t ny, std::size_t nz, double dx, double dy,
                                           double dz, double tolerance, std::size_t maxIterations) {
    const std::size_t n = nx * ny * nz;
    std::vector<double> p = initialGuess;
    const std::vector<double> ap = cpuReferencePoisson(p, nx, ny, nz, dx, dy, dz);
    std::vector<double> residual(n);
    for (std::size_t idx = 0; idx < n; ++idx) {
        residual[idx] = (idx == 0) ? 0.0 : rhs[idx] - ap[idx];
    }
    std::vector<double> direction = residual;
    double residualDotResidual = dotProduct(residual, residual);

    std::size_t iteration = 0;
    bool brokeDown = false;
    for (; iteration < maxIterations; ++iteration) {
        if (std::sqrt(residualDotResidual) < tolerance) {
            break;
        }
        const std::vector<double> ad = cpuReferencePoisson(direction, nx, ny, nz, dx, dy, dz);
        const double directionDotAd = dotProduct(direction, ad);
        if (directionDotAd == 0.0) {
            brokeDown = true;
            break;
        }
        const double alpha = residualDotResidual / directionDotAd;
        for (std::size_t idx = 1; idx < n; ++idx) {
            p[idx] += alpha * direction[idx];
        }
        for (std::size_t idx = 0; idx < n; ++idx) {
            residual[idx] = (idx == 0) ? 0.0 : residual[idx] - alpha * ad[idx];
        }
        const double newResidualDotResidual = dotProduct(residual, residual);
        const double beta = newResidualDotResidual / residualDotResidual;
        for (std::size_t idx = 0; idx < n; ++idx) {
            direction[idx] = (idx == 0) ? 0.0 : residual[idx] + beta * direction[idx];
        }
        residualDotResidual = newResidualDotResidual;
    }

    CpuCgResult result;
    result.pressure = std::move(p);
    result.iterations = iteration;
    result.residualNorm = std::sqrt(residualDotResidual);
    result.brokeDown = brokeDown;
    return result;
}

// Check A: the GPU's own returned pressure field independently
// re-satisfies residual = rhs - A*p when that residual is recomputed from
// scratch on the CPU (a fresh stencil apply plus a plain sequential dot,
// not the GPU's own incrementally-maintained residual) -- this validates
// "is this actually a solution to the stated problem", independent of
// whether the GPU's internal iteration trajectory matched a CPU one bit
// for bit, which the class's own header explains is not expected (a
// parallel dot-product reduction sums in a different, though fixed and
// reproducible, order than a sequential one).
void testConjugateGradientSolverCudaResidualIsSmallIndependently() {
    const std::size_t nx = 12;
    const std::size_t ny = 10;
    const std::size_t nz = 8;
    const double dx = 0.1;
    const double dy = 0.15;
    const double dz = 0.2;

    ConjugateGradientSolverCuda solver(nx, ny, nz, dx, dy, dz);
    if (!solver.available()) {
        std::printf("  [aether_gpu_tests] nenhum dispositivo CUDA disponivel -- pulando "
                     "ConjugateGradientSolverCuda.\n");
        return;
    }

    const std::vector<double> rhs = testField(nx, ny, nz);
    std::vector<double> initialGuess(nx * ny * nz);
    for (std::size_t i = 0; i < initialGuess.size(); ++i) {
        initialGuess[i] = 1.0 + 0.3 * std::sin(0.4 * static_cast<double>(i));
    }

    const auto result = solver.solve(rhs, initialGuess);
    AETHER_CHECK(!result.brokeDown);
    AETHER_CHECK(result.converged);

    const std::vector<double> ap = cpuReferencePoisson(result.pressure, nx, ny, nz, dx, dy, dz);
    double residualNormSquared = 0.0;
    for (std::size_t idx = 0; idx < ap.size(); ++idx) {
        const double r = (idx == 0) ? 0.0 : rhs[idx] - ap[idx];
        residualNormSquared += r * r;
    }
    const double independentResidualNorm = std::sqrt(residualNormSquared);
    std::printf("  [aether_gpu_tests] CG residente: %zu iteracoes, residual reportado pela GPU=%.3e, "
                "residual reverificado na CPU=%.3e\n",
                result.iterations, result.residualNorm, independentResidualNorm);
    // Medido, nao suposto: 9.953e-11 nesta grade -- a re-soma sequencial na
    // CPU do residuo (a partir da pressao final da GPU) fica na mesma
    // ordem de grandeza da tolerancia de convergencia (1e-10) que a
    // propria GPU usou internamente, com folga de ~10x acima do valor
    // medido.
    AETHER_CHECK(independentResidualNorm < 1e-9);
}

// Check B: the pinned cell (index 0) is never touched by
// updatePressureKernel, so Result.pressure[0] must equal
// initialGuess[0] bit-for-bit regardless of tolerance/iterations/rhs --
// cheap, always-applicable, and catches immediately if any one of the
// three update kernels ever loses its own idx==0 guard. initialGuess[0]
// is deliberately a distinctive, non-zero value here so the check isn't
// vacuously true.
void testConjugateGradientSolverCudaPinnedCellMatchesInitialGuess() {
    const std::size_t nx = 5;
    const std::size_t ny = 5;
    const std::size_t nz = 5;
    ConjugateGradientSolverCuda solver(nx, ny, nz, 1.0, 1.0, 1.0);
    if (!solver.available()) {
        return; // already reported above
    }
    const std::vector<double> rhs = testField(nx, ny, nz);
    std::vector<double> initialGuess(nx * ny * nz, 2.0);
    initialGuess[0] = -7.25;

    const auto result = solver.solve(rhs, initialGuess);
    AETHER_CHECK(result.pressure[0] == initialGuess[0]);
}

// Check C: GPU trajectory vs. an independent CPU reference CG loop on the
// same problem, tolerance measured and printed rather than assumed. The
// printed diagnostic is what determines the AETHER_CHECK threshold below
// -- not the other way around.
void testConjugateGradientSolverCudaTrajectoryMatchesCpuReference() {
    const std::size_t nx = 10;
    const std::size_t ny = 9;
    const std::size_t nz = 8;
    const double dx = 0.12;
    const double dy = 0.09;
    const double dz = 0.15;

    ConjugateGradientSolverCuda solver(nx, ny, nz, dx, dy, dz);
    if (!solver.available()) {
        return; // already reported above
    }

    const std::vector<double> rhs = testField(nx, ny, nz);
    const std::vector<double> initialGuess(nx * ny * nz, 0.0);

    const auto gpuResult = solver.solve(rhs, initialGuess);
    const auto cpuResult =
        cpuReferenceConjugateGradient(rhs, initialGuess, nx, ny, nz, dx, dy, dz, 1e-10, nx * ny * nz);

    double maxDiff = 0.0;
    for (std::size_t i = 0; i < gpuResult.pressure.size(); ++i) {
        maxDiff = std::max(maxDiff, std::fabs(gpuResult.pressure[i] - cpuResult.pressure[i]));
    }
    std::printf("  [aether_gpu_tests] CG: GPU %zu iteracoes (residual %.3e) vs CPU %zu iteracoes "
                "(residual %.3e), maxDiff no campo final=%.3e\n",
                gpuResult.iterations, gpuResult.residualNorm, cpuResult.iterations, cpuResult.residualNorm,
                maxDiff);

    AETHER_CHECK(gpuResult.converged);
    AETHER_CHECK(!gpuResult.brokeDown);
    AETHER_CHECK(!cpuResult.brokeDown);
    // Medido, nao suposto: 117 iteracoes dos dois lados (identico) e
    // maxDiff=2.398e-14 nesta grade -- a diferenca de ordem de soma do
    // produto interno paralelo acaba sendo desprezivel aqui, bem abaixo
    // do epsilon de double (~2.2e-16 por operacao, acumulado por poucas
    // dezenas de iteracoes). Limiar com ~4 ordens de grandeza de folga.
    AETHER_CHECK(maxDiff < 1e-10);
}

// n=1: the CPU loop's own structure shows residual[0] is forced to 0
// regardless of grid size, so residualDotResidual is exactly 0 before any
// per-iteration kernel ever runs on `direction` -- the loop must exit at
// iteration 0, having only computed the initial Ap = applyLaplacian(p).
void testConjugateGradientSolverCudaHandlesSingleCellGrid() {
    ConjugateGradientSolverCuda solver(1, 1, 1, 1.0, 1.0, 1.0);
    if (!solver.available()) {
        return; // already reported above
    }
    const std::vector<double> rhs = {5.0};
    const std::vector<double> initialGuess = {3.0};

    const auto result = solver.solve(rhs, initialGuess);
    AETHER_CHECK(result.iterations == 0);
    AETHER_CHECK(result.converged);
    AETHER_CHECK(!result.brokeDown);
    AETHER_CHECK(result.pressure.size() == 1);
    AETHER_CHECK(result.pressure[0] == 3.0);
}

// A grid whose cell count (120) is not a multiple of the 256-thread
// reduction/AXPY block size -- exercises the grid-stride bound in
// dotProductPartialKernel/sumReduceKernel and the 1D bounds check in the
// update kernels simultaneously. Same grid shape
// testPoissonOperatorCudaMatchesCpuExactly() already uses above.
void testConjugateGradientSolverCudaHandlesNonMultipleOfBlockSizeGrid() {
    const std::size_t nx = 6;
    const std::size_t ny = 5;
    const std::size_t nz = 4;
    const double dx = 0.1;
    const double dy = 0.2;
    const double dz = 0.05;

    ConjugateGradientSolverCuda solver(nx, ny, nz, dx, dy, dz);
    if (!solver.available()) {
        return; // already reported above
    }
    const std::vector<double> rhs = testField(nx, ny, nz);
    const std::vector<double> initialGuess(nx * ny * nz, 0.0);

    const auto result = solver.solve(rhs, initialGuess);
    AETHER_CHECK(result.converged);
    AETHER_CHECK(!result.brokeDown);
    AETHER_CHECK(result.pressure[0] == initialGuess[0]);
    std::printf("  [aether_gpu_tests] CG grade 6x5x4 (nao multiplo do bloco de 256): %zu iteracoes, "
                "residual=%.3e\n",
                result.iterations, result.residualNorm);
}

// --- MomentumPredictorCuda: independent CPU reference, Central scheme
// only (see MomentumPredictorCuda.hpp for why that's the only scheme that
// needs porting), reduced from StaggeredCavityBase3D::computeMomentumPredictor
// the same way the GPU kernels themselves were -- transcribed here rather
// than calling into aether_solver (same discipline as cpuReferencePoisson
// above; this test target links only aether_gpu + aether_testing).

std::size_t refIndexU(std::size_t i, std::size_t j, std::size_t k, std::size_t nx, std::size_t ny) {
    return i + j * (nx + 1) + k * (nx + 1) * ny;
}
std::size_t refIndexV(std::size_t i, std::size_t j, std::size_t k, std::size_t nx, std::size_t ny) {
    return i + j * nx + k * nx * (ny + 1);
}
std::size_t refIndexW(std::size_t i, std::size_t j, std::size_t k, std::size_t nx, std::size_t ny) {
    return i + j * nx + k * nx * ny;
}
std::size_t refIndexP(std::size_t i, std::size_t j, std::size_t k, std::size_t nx, std::size_t ny) {
    return i + j * nx + k * nx * ny;
}

double refLidVelocityAt(double x, double y, double lengthX, double lengthY, double lidVelocity) {
    constexpr double kPi = 3.14159265358979323846;
    const double sx = std::sin(kPi * x / lengthX);
    const double sy = std::sin(kPi * y / lengthY);
    return lidVelocity * sx * sx * sy * sy;
}

double refUAt(const std::vector<double>& u, long long i, long long j, long long k, std::size_t nx,
              std::size_t ny, std::size_t nz, double dx, double dy, double lengthX, double lengthY,
              double lidVelocity) {
    const auto ii = static_cast<std::size_t>(i);
    if (j < 0) {
        return -u[refIndexU(ii, 0, static_cast<std::size_t>(k), nx, ny)];
    }
    if (j >= static_cast<long long>(ny)) {
        return -u[refIndexU(ii, ny - 1, static_cast<std::size_t>(k), nx, ny)];
    }
    if (k < 0) {
        return -u[refIndexU(ii, static_cast<std::size_t>(j), 0, nx, ny)];
    }
    if (k >= static_cast<long long>(nz)) {
        const double x = static_cast<double>(i) * dx;
        const double y = (static_cast<double>(j) + 0.5) * dy;
        return 2.0 * refLidVelocityAt(x, y, lengthX, lengthY, lidVelocity) -
               u[refIndexU(ii, static_cast<std::size_t>(j), nz - 1, nx, ny)];
    }
    return u[refIndexU(ii, static_cast<std::size_t>(j), static_cast<std::size_t>(k), nx, ny)];
}

double refVAt(const std::vector<double>& v, long long i, long long j, long long k, std::size_t nx,
              std::size_t ny, std::size_t nz) {
    const auto jj = static_cast<std::size_t>(j);
    if (i < 0) {
        return -v[refIndexV(0, jj, static_cast<std::size_t>(k), nx, ny)];
    }
    if (i >= static_cast<long long>(nx)) {
        return -v[refIndexV(nx - 1, jj, static_cast<std::size_t>(k), nx, ny)];
    }
    if (k < 0) {
        return -v[refIndexV(static_cast<std::size_t>(i), jj, 0, nx, ny)];
    }
    if (k >= static_cast<long long>(nz)) {
        return -v[refIndexV(static_cast<std::size_t>(i), jj, nz - 1, nx, ny)];
    }
    return v[refIndexV(static_cast<std::size_t>(i), jj, static_cast<std::size_t>(k), nx, ny)];
}

double refWAt(const std::vector<double>& w, long long i, long long j, long long k, std::size_t nx,
              std::size_t ny, std::size_t nz) {
    const auto kk = static_cast<std::size_t>(k);
    if (i < 0) {
        return -w[refIndexW(0, static_cast<std::size_t>(j), kk, nx, ny)];
    }
    if (i >= static_cast<long long>(nx)) {
        return -w[refIndexW(nx - 1, static_cast<std::size_t>(j), kk, nx, ny)];
    }
    if (j < 0) {
        return -w[refIndexW(static_cast<std::size_t>(i), 0, kk, nx, ny)];
    }
    if (j >= static_cast<long long>(ny)) {
        return -w[refIndexW(static_cast<std::size_t>(i), ny - 1, kk, nx, ny)];
    }
    return w[refIndexW(static_cast<std::size_t>(i), static_cast<std::size_t>(j), kk, nx, ny)];
}

double refNutAt(const std::vector<double>& nut, long long i, long long j, long long k, std::size_t nx,
                std::size_t ny, std::size_t nz) {
    if (nut.empty()) {
        return 0.0;
    }
    if (i < 0 || i >= static_cast<long long>(nx) || j < 0 || j >= static_cast<long long>(ny) || k < 0 ||
        k >= static_cast<long long>(nz)) {
        return 0.0;
    }
    return nut[refIndexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k),
                          nx, ny)];
}

struct MomentumPredictorResult {
    std::vector<double> uStar;
    std::vector<double> vStar;
    std::vector<double> wStar;
};

MomentumPredictorResult cpuReferenceMomentumPredictor(const std::vector<double>& u, const std::vector<double>& v,
                                                       const std::vector<double>& w,
                                                       const std::vector<double>& nut, std::size_t nx,
                                                       std::size_t ny, std::size_t nz, double dx, double dy,
                                                       double dz, double lengthX, double lengthY,
                                                       double viscosity, double lidVelocity, double dt) {
    MomentumPredictorResult result;
    result.uStar = u;
    result.vStar = v;
    result.wStar = w;

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 1; i < nx; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double uHere = refUAt(u, li, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
                const double uIp1 = refUAt(u, li + 1, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
                const double uIm1 = refUAt(u, li - 1, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
                const double uCenterI = 0.5 * (uHere + uIp1);
                const double uCenterIm1 = 0.5 * (uIm1 + uHere);
                const double duudx = (uCenterI * uCenterI - uCenterIm1 * uCenterIm1) / dx;

                const double uJp1 = refUAt(u, li, lj + 1, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
                const double uJm1 = refUAt(u, li, lj - 1, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
                const double uEdgeJp = 0.5 * (uHere + uJp1);
                const double uEdgeJm = 0.5 * (uJm1 + uHere);
                const double vEdgeIJp =
                    0.5 * (refVAt(v, li - 1, lj + 1, lk, nx, ny, nz) + refVAt(v, li, lj + 1, lk, nx, ny, nz));
                const double vEdgeIJ = 0.5 * (refVAt(v, li - 1, lj, lk, nx, ny, nz) + refVAt(v, li, lj, lk, nx, ny, nz));
                const double duvdy = (vEdgeIJp * uEdgeJp - vEdgeIJ * uEdgeJm) / dy;

                const double uKp1 = refUAt(u, li, lj, lk + 1, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
                const double uKm1 = refUAt(u, li, lj, lk - 1, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
                const double uEdgeKp = 0.5 * (uHere + uKp1);
                const double uEdgeKm = 0.5 * (uKm1 + uHere);
                const double wEdgeIKp =
                    0.5 * (refWAt(w, li - 1, lj, lk + 1, nx, ny, nz) + refWAt(w, li, lj, lk + 1, nx, ny, nz));
                const double wEdgeIK = 0.5 * (refWAt(w, li - 1, lj, lk, nx, ny, nz) + refWAt(w, li, lj, lk, nx, ny, nz));
                const double duwdz = (wEdgeIKp * uEdgeKp - wEdgeIK * uEdgeKm) / dz;

                const double gammaE = viscosity + refNutAt(nut, li, lj, lk, nx, ny, nz);
                const double gammaW = viscosity + refNutAt(nut, li - 1, lj, lk, nx, ny, nz);
                const double gammaT = viscosity + 0.5 * (refNutAt(nut, li - 1, lj, lk, nx, ny, nz) +
                                                          refNutAt(nut, li, lj, lk, nx, ny, nz));
                const double diffusionU = (gammaE * (uIp1 - uHere) - gammaW * (uHere - uIm1)) / (dx * dx) +
                                           gammaT * (uJp1 - 2.0 * uHere + uJm1) / (dy * dy) +
                                           gammaT * (uKp1 - 2.0 * uHere + uKm1) / (dz * dz);

                result.uStar[refIndexU(i, j, k, nx, ny)] = uHere + dt * (-(duudx + duvdy + duwdz) + diffusionU);
            }
        }
    }

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 1; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double vHere = refVAt(v, li, lj, lk, nx, ny, nz);
                const double vJp1 = refVAt(v, li, lj + 1, lk, nx, ny, nz);
                const double vJm1 = refVAt(v, li, lj - 1, lk, nx, ny, nz);
                const double vCenterJ = 0.5 * (vHere + vJp1);
                const double vCenterJm1 = 0.5 * (vJm1 + vHere);
                const double dvvdy = (vCenterJ * vCenterJ - vCenterJm1 * vCenterJm1) / dy;

                const double vKp1 = refVAt(v, li, lj, lk + 1, nx, ny, nz);
                const double vKm1 = refVAt(v, li, lj, lk - 1, nx, ny, nz);
                const double vEdgeKp = 0.5 * (vHere + vKp1);
                const double vEdgeKm = 0.5 * (vKm1 + vHere);
                const double wEdgeJKp =
                    0.5 * (refWAt(w, li, lj - 1, lk + 1, nx, ny, nz) + refWAt(w, li, lj, lk + 1, nx, ny, nz));
                const double wEdgeJK = 0.5 * (refWAt(w, li, lj - 1, lk, nx, ny, nz) + refWAt(w, li, lj, lk, nx, ny, nz));
                const double dvwdz = (wEdgeJKp * vEdgeKp - wEdgeJK * vEdgeKm) / dz;

                const double vIp1 = refVAt(v, li + 1, lj, lk, nx, ny, nz);
                const double vIm1 = refVAt(v, li - 1, lj, lk, nx, ny, nz);
                const double vEdgeIp = 0.5 * (vHere + vIp1);
                const double vEdgeIm = 0.5 * (vIm1 + vHere);
                const double uEdgeJIp = 0.5 * (refUAt(u, li + 1, lj - 1, lk, nx, ny, nz, dx, dy, lengthX, lengthY,
                                                       lidVelocity) +
                                                refUAt(u, li + 1, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY,
                                                       lidVelocity));
                const double uEdgeJI = 0.5 * (refUAt(u, li, lj - 1, lk, nx, ny, nz, dx, dy, lengthX, lengthY,
                                                      lidVelocity) +
                                               refUAt(u, li, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY,
                                                      lidVelocity));
                const double dvudx = (uEdgeJIp * vEdgeIp - uEdgeJI * vEdgeIm) / dx;

                const double gammaN = viscosity + refNutAt(nut, li, lj, lk, nx, ny, nz);
                const double gammaS = viscosity + refNutAt(nut, li, lj - 1, lk, nx, ny, nz);
                const double gammaT = viscosity + 0.5 * (refNutAt(nut, li, lj - 1, lk, nx, ny, nz) +
                                                          refNutAt(nut, li, lj, lk, nx, ny, nz));
                const double diffusionV = (gammaN * (vJp1 - vHere) - gammaS * (vHere - vJm1)) / (dy * dy) +
                                           gammaT * (vIp1 - 2.0 * vHere + vIm1) / (dx * dx) +
                                           gammaT * (vKp1 - 2.0 * vHere + vKm1) / (dz * dz);

                result.vStar[refIndexV(i, j, k, nx, ny)] = vHere + dt * (-(dvvdy + dvwdz + dvudx) + diffusionV);
            }
        }
    }

    for (std::size_t k = 1; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double wHere = refWAt(w, li, lj, lk, nx, ny, nz);
                const double wKp1 = refWAt(w, li, lj, lk + 1, nx, ny, nz);
                const double wKm1 = refWAt(w, li, lj, lk - 1, nx, ny, nz);
                const double wCenterK = 0.5 * (wHere + wKp1);
                const double wCenterKm1 = 0.5 * (wKm1 + wHere);
                const double dwwdz = (wCenterK * wCenterK - wCenterKm1 * wCenterKm1) / dz;

                const double wIp1 = refWAt(w, li + 1, lj, lk, nx, ny, nz);
                const double wIm1 = refWAt(w, li - 1, lj, lk, nx, ny, nz);
                const double wEdgeIp = 0.5 * (wHere + wIp1);
                const double wEdgeIm = 0.5 * (wIm1 + wHere);
                const double uEdgeKIp = 0.5 * (refUAt(u, li + 1, lj, lk - 1, nx, ny, nz, dx, dy, lengthX, lengthY,
                                                       lidVelocity) +
                                                refUAt(u, li + 1, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY,
                                                       lidVelocity));
                const double uEdgeKI = 0.5 * (refUAt(u, li, lj, lk - 1, nx, ny, nz, dx, dy, lengthX, lengthY,
                                                      lidVelocity) +
                                               refUAt(u, li, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY,
                                                      lidVelocity));
                const double dwudx = (uEdgeKIp * wEdgeIp - uEdgeKI * wEdgeIm) / dx;

                const double wJp1 = refWAt(w, li, lj + 1, lk, nx, ny, nz);
                const double wJm1 = refWAt(w, li, lj - 1, lk, nx, ny, nz);
                const double wEdgeJp = 0.5 * (wHere + wJp1);
                const double wEdgeJm = 0.5 * (wJm1 + wHere);
                const double vEdgeKJp =
                    0.5 * (refVAt(v, li, lj + 1, lk - 1, nx, ny, nz) + refVAt(v, li, lj + 1, lk, nx, ny, nz));
                const double vEdgeKJ = 0.5 * (refVAt(v, li, lj, lk - 1, nx, ny, nz) + refVAt(v, li, lj, lk, nx, ny, nz));
                const double dwvdy = (vEdgeKJp * wEdgeJp - vEdgeKJ * wEdgeJm) / dy;

                const double gammaF = viscosity + refNutAt(nut, li, lj, lk, nx, ny, nz);
                const double gammaB = viscosity + refNutAt(nut, li, lj, lk - 1, nx, ny, nz);
                const double gammaT = viscosity + 0.5 * (refNutAt(nut, li, lj, lk - 1, nx, ny, nz) +
                                                          refNutAt(nut, li, lj, lk, nx, ny, nz));
                const double diffusionW = (gammaF * (wKp1 - wHere) - gammaB * (wHere - wKm1)) / (dz * dz) +
                                           gammaT * (wIp1 - 2.0 * wHere + wIm1) / (dx * dx) +
                                           gammaT * (wJp1 - 2.0 * wHere + wJm1) / (dy * dy);

                result.wStar[refIndexW(i, j, k, nx, ny)] = wHere + dt * (-(dwwdz + dwudx + dwvdy) + diffusionW);
            }
        }
    }

    return result;
}

// Deterministic field with its own distinct coefficients -- deliberately
// NOT the same generator relabeled per component. A u<->v<->w mixup
// anywhere in the port (e.g. a cross term reading nutAt(i,j-1,k) where
// the w-kernel needs nutAt(i,j,k-1)) must produce a guaranteed-different
// result, not a coincidentally-plausible one from an equally-shaped
// generator.
std::vector<double> generateTestField(std::size_t size, double c1, double c2, double c3, double c4) {
    std::vector<double> field(size);
    for (std::size_t idx = 0; idx < size; ++idx) {
        const double x = static_cast<double>(idx);
        field[idx] = c1 * std::sin(c2 * x) + c3 * std::cos(c4 * x) + 0.01 * x;
    }
    return field;
}

void checkMomentumResultsMatchExactly(const MomentumPredictorCuda::Result& gpuResult,
                                       const MomentumPredictorResult& cpuResult, const char* label) {
    AETHER_CHECK(gpuResult.uStar.size() == cpuResult.uStar.size());
    AETHER_CHECK(gpuResult.vStar.size() == cpuResult.vStar.size());
    AETHER_CHECK(gpuResult.wStar.size() == cpuResult.wStar.size());

    double maxDiff = 0.0;
    std::size_t mismatches = 0;
    const auto compare = [&](const std::vector<double>& a, const std::vector<double>& b) {
        for (std::size_t i = 0; i < a.size(); ++i) {
            maxDiff = std::max(maxDiff, std::fabs(a[i] - b[i]));
            if (a[i] != b[i]) {
                ++mismatches;
            }
        }
    };
    compare(gpuResult.uStar, cpuResult.uStar);
    compare(gpuResult.vStar, cpuResult.vStar);
    compare(gpuResult.wStar, cpuResult.wStar);
    std::printf("  [aether_gpu_tests] MomentumPredictorCuda (%s): maxDiff=%.3e, mismatches=%zu\n", label, maxDiff,
                mismatches);

    for (std::size_t i = 0; i < gpuResult.uStar.size(); ++i) {
        AETHER_CHECK(gpuResult.uStar[i] == cpuResult.uStar[i]);
    }
    for (std::size_t i = 0; i < gpuResult.vStar.size(); ++i) {
        AETHER_CHECK(gpuResult.vStar[i] == cpuResult.vStar[i]);
    }
    for (std::size_t i = 0; i < gpuResult.wStar.size(); ++i) {
        AETHER_CHECK(gpuResult.wStar[i] == cpuResult.wStar[i]);
    }
}

// Non-multiple-of-block-size grid (6x5x4, same as the Poisson/CG tests
// above), laminar (nut empty).
void testMomentumPredictorCudaMatchesCpuExactlyLaminar() {
    const std::size_t nx = 6, ny = 5, nz = 4;
    const double lengthX = 1.3, lengthY = 0.9, lengthZ = 1.7;
    const double viscosity = 0.05, lidVelocity = 1.2, dt = 0.01;
    const double dx = lengthX / static_cast<double>(nx);
    const double dy = lengthY / static_cast<double>(ny);
    const double dz = lengthZ / static_cast<double>(nz);

    MomentumPredictorCuda predictor(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity);
    if (!predictor.available()) {
        std::printf("  [aether_gpu_tests] nenhum dispositivo CUDA disponivel -- pulando "
                     "MomentumPredictorCuda.\n");
        return;
    }

    const std::vector<double> u = generateTestField((nx + 1) * ny * nz, 0.70, 0.31, 0.40, 0.21);
    const std::vector<double> v = generateTestField(nx * (ny + 1) * nz, 0.50, 0.19, 0.30, 0.37);
    const std::vector<double> w = generateTestField(nx * ny * (nz + 1), 0.60, 0.27, 0.20, 0.13);
    const std::vector<double> nut; // laminar

    const auto gpuResult = predictor.predict(u, v, w, nut, dt);
    const auto cpuResult =
        cpuReferenceMomentumPredictor(u, v, w, nut, nx, ny, nz, dx, dy, dz, lengthX, lengthY, viscosity,
                                       lidVelocity, dt);
    checkMomentumResultsMatchExactly(gpuResult, cpuResult, "laminar, grade 6x5x4");
}

// Same grid, non-zero deterministic nut -- confirms gammaE/gammaW/
// gammaTransverse port correctly, not just the laminar (nut==0) path.
void testMomentumPredictorCudaMatchesCpuExactlyTurbulent() {
    const std::size_t nx = 6, ny = 5, nz = 4;
    const double lengthX = 1.3, lengthY = 0.9, lengthZ = 1.7;
    const double viscosity = 0.05, lidVelocity = 1.2, dt = 0.01;
    const double dx = lengthX / static_cast<double>(nx);
    const double dy = lengthY / static_cast<double>(ny);
    const double dz = lengthZ / static_cast<double>(nz);

    MomentumPredictorCuda predictor(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity);
    if (!predictor.available()) {
        return; // already reported above
    }

    const std::vector<double> u = generateTestField((nx + 1) * ny * nz, 0.70, 0.31, 0.40, 0.21);
    const std::vector<double> v = generateTestField(nx * (ny + 1) * nz, 0.50, 0.19, 0.30, 0.37);
    const std::vector<double> w = generateTestField(nx * ny * (nz + 1), 0.60, 0.27, 0.20, 0.13);
    const std::vector<double> nut = generateTestField(nx * ny * nz, 0.03, 0.15, 0.02, 0.22);

    const auto gpuResult = predictor.predict(u, v, w, nut, dt);
    const auto cpuResult =
        cpuReferenceMomentumPredictor(u, v, w, nut, nx, ny, nz, dx, dy, dz, lengthX, lengthY, viscosity,
                                       lidVelocity, dt);
    checkMomentumResultsMatchExactly(gpuResult, cpuResult, "turbulento, grade 6x5x4");
}

// Black-box invariant: every boundary face (i in {0,nx} for u, j in
// {0,ny} for v, k in {0,nz} for w) must equal the input bit-for-bit,
// independent of whether the interior formula is even correct -- catches
// an off-by-one in a kernel's own boundary guard in isolation.
void testMomentumPredictorCudaPreservesBoundaryFaces() {
    const std::size_t nx = 5, ny = 5, nz = 5;
    const double lengthX = 1.0, lengthY = 1.0, lengthZ = 1.0;
    const double viscosity = 0.1, lidVelocity = 1.0, dt = 0.01;

    MomentumPredictorCuda predictor(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity);
    if (!predictor.available()) {
        return; // already reported above
    }

    const std::vector<double> u = generateTestField((nx + 1) * ny * nz, 0.90, 0.44, 0.12, 0.33);
    const std::vector<double> v = generateTestField(nx * (ny + 1) * nz, 0.80, 0.29, 0.22, 0.18);
    const std::vector<double> w = generateTestField(nx * ny * (nz + 1), 0.70, 0.17, 0.31, 0.26);
    const std::vector<double> nut;

    const auto result = predictor.predict(u, v, w, nut, dt);

    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            AETHER_CHECK(result.uStar[refIndexU(0, j, k, nx, ny)] == u[refIndexU(0, j, k, nx, ny)]);
            AETHER_CHECK(result.uStar[refIndexU(nx, j, k, nx, ny)] == u[refIndexU(nx, j, k, nx, ny)]);
        }
    }
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t i = 0; i < nx; ++i) {
            AETHER_CHECK(result.vStar[refIndexV(i, 0, k, nx, ny)] == v[refIndexV(i, 0, k, nx, ny)]);
            AETHER_CHECK(result.vStar[refIndexV(i, ny, k, nx, ny)] == v[refIndexV(i, ny, k, nx, ny)]);
        }
    }
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            AETHER_CHECK(result.wStar[refIndexW(i, j, 0, nx, ny)] == w[refIndexW(i, j, 0, nx, ny)]);
            AETHER_CHECK(result.wStar[refIndexW(i, j, nz, nx, ny)] == w[refIndexW(i, j, nz, nx, ny)]);
        }
    }
}

// The nut device buffer is reused across calls on the same instance --
// a turbulent call followed by a laminar call on the SAME predictor must
// not let the turbulent field leak into the laminar result (see
// MomentumPredictorCuda::predict()'s own comment on why the laminar path
// explicitly re-zeros the buffer instead of relying on it starting at
// zero).
void testMomentumPredictorCudaReZerosNutBetweenCalls() {
    const std::size_t nx = 5, ny = 5, nz = 5;
    const double lengthX = 1.0, lengthY = 1.0, lengthZ = 1.0;
    const double viscosity = 0.1, lidVelocity = 1.0, dt = 0.01;
    const double dx = lengthX / static_cast<double>(nx);
    const double dy = lengthY / static_cast<double>(ny);
    const double dz = lengthZ / static_cast<double>(nz);

    MomentumPredictorCuda predictor(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity);
    if (!predictor.available()) {
        return; // already reported above
    }

    const std::vector<double> u = generateTestField((nx + 1) * ny * nz, 0.65, 0.34, 0.11, 0.24);
    const std::vector<double> v = generateTestField(nx * (ny + 1) * nz, 0.55, 0.21, 0.19, 0.28);
    const std::vector<double> w = generateTestField(nx * ny * (nz + 1), 0.45, 0.16, 0.23, 0.31);
    const std::vector<double> nut = generateTestField(nx * ny * nz, 0.03, 0.12, 0.02, 0.19);

    (void)predictor.predict(u, v, w, nut, dt); // turbulent call first, dirties the reused nut buffer

    const std::vector<double> empty;
    const auto gpuLaminar = predictor.predict(u, v, w, empty, dt); // laminar call second, same instance
    const auto cpuLaminar =
        cpuReferenceMomentumPredictor(u, v, w, empty, nx, ny, nz, dx, dy, dz, lengthX, lengthY, viscosity,
                                       lidVelocity, dt);
    checkMomentumResultsMatchExactly(gpuLaminar, cpuLaminar, "laminar apos turbulento, mesma instancia");
}

} // namespace

int main() {
    testPoissonOperatorCudaMatchesCpuExactly();
    testPoissonOperatorCudaPinnedCellIsIdentity();
    testConjugateGradientSolverCudaResidualIsSmallIndependently();
    testConjugateGradientSolverCudaPinnedCellMatchesInitialGuess();
    testConjugateGradientSolverCudaTrajectoryMatchesCpuReference();
    testConjugateGradientSolverCudaHandlesSingleCellGrid();
    testConjugateGradientSolverCudaHandlesNonMultipleOfBlockSizeGrid();
    testMomentumPredictorCudaMatchesCpuExactlyLaminar();
    testMomentumPredictorCudaMatchesCpuExactlyTurbulent();
    testMomentumPredictorCudaPreservesBoundaryFaces();
    testMomentumPredictorCudaReZerosNutBetweenCalls();
    std::printf("aether_gpu_tests: OK\n");
    return 0;
}
