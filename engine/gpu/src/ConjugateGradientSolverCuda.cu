#include "aether/gpu/ConjugateGradientSolverCuda.hpp"

#include "detail/PoissonStencilCuda.cuh"

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace aether::gpu {

namespace {

constexpr int kBlockSize1D = 256;
constexpr std::size_t kMaxReductionBlocks = 4096;

// The stencil apply, 3D launch -- same (8,8,8) block shape as
// PoissonOperatorCuda's own poissonKernel, same shared per-cell formula
// (detail/PoissonStencilCuda.cuh). Orchestration (launch config, which
// buffer feeds which) stays independent per file; only the formula is
// shared -- see that header's own comment for why.
__global__ void applyLaplacianKernel(const double* x, double* result, std::size_t nx, std::size_t ny,
                                      std::size_t nz, double ax, double ay, double az, double weightTotal) {
    const auto i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto j = static_cast<std::size_t>(blockIdx.y) * blockDim.y + threadIdx.y;
    const auto k = static_cast<std::size_t>(blockIdx.z) * blockDim.z + threadIdx.z;
    if (i >= nx || j >= ny || k >= nz) {
        return;
    }
    const std::size_t idx = i + j * nx + k * nx * ny;
    result[idx] = detail::poissonStencilValue(x, i, j, k, nx, ny, nz, ax, ay, az, weightTotal);
}

// Sums the 32 lanes of one warp via shuffle -- no shared memory needed at
// this level. Valid in every lane after the loop (a butterfly-style
// shuffle-down reduction), but only lane 0's value is used by the caller.
__device__ __forceinline__ double warpReduceSum(double v) {
    for (int offset = 16; offset > 0; offset >>= 1) {
        v += __shfl_down_sync(0xffffffffu, v, offset);
    }
    return v;
}

// Reduces one block's threads to a single sum, valid only in thread 0.
// Two-level: each warp reduces itself via warpReduceSum, the (at most 32,
// for blockDim<=1024) per-warp sums are written to shared memory, then
// warp 0 reduces those the same way. Deliberately not a single-pass
// atomicAdd-based reduction -- see dotProductPartialKernel's own comment
// for why determinism here matters to this class's own validation story.
__device__ double blockReduceSum(double v) {
    __shared__ double warpSums[32]; // up to 1024 threads/block -> 32 warps
    const int lane = threadIdx.x & 31;
    const int warpId = threadIdx.x >> 5;

    v = warpReduceSum(v);
    if (lane == 0) {
        warpSums[warpId] = v;
    }
    __syncthreads();

    const int numWarps = (blockDim.x + 31) / 32;
    v = (threadIdx.x < numWarps) ? warpSums[threadIdx.x] : 0.0;
    if (warpId == 0) {
        v = warpReduceSum(v);
    }
    return v;
}

// First pass of the dot product: each block grid-strides over the full
// array and writes exactly one partial sum. Grid-stride (not "one thread
// covers a fixed slice") makes this correct for any n regardless of how
// numBlocks was chosen, including n not a multiple of the block/grid
// size.
//
// **Why two passes with a fixed number of blocks, not one pass with
// atomicAdd**: an atomicAdd-based reduction sums whichever block finishes
// first into a running total, in scheduling-dependent order that varies
// run to run on identical inputs -- non-deterministic. This class's own
// validation deliberately measures the *legitimate* difference between a
// parallel reduction's summation order and the CPU's sequential one (see
// the class header); a non-deterministic GPU answer would make that
// measurement a moving target instead of a fixed, explainable number.
// This design keeps the summation order fixed by block index for a given
// launch configuration, so the GPU's own answer is reproducible run to
// run, and only genuinely differs from the CPU's for the one reason
// that's expected to.
__global__ void dotProductPartialKernel(const double* a, const double* b, double* partialSums, std::size_t n) {
    const std::size_t start = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const std::size_t stride = static_cast<std::size_t>(gridDim.x) * blockDim.x;
    double sum = 0.0;
    for (std::size_t i = start; i < n; i += stride) {
        sum += a[i] * b[i];
    }
    sum = blockReduceSum(sum);
    if (threadIdx.x == 0) {
        partialSums[blockIdx.x] = sum;
    }
}

// Second pass: a single block sums the (at most kMaxReductionBlocks, so
// always small) partial sums from the first pass into one scalar.
__global__ void sumReduceKernel(const double* values, std::size_t count, double* result) {
    double sum = 0.0;
    for (std::size_t i = threadIdx.x; i < count; i += blockDim.x) {
        sum += values[i];
    }
    sum = blockReduceSum(sum);
    if (threadIdx.x == 0) {
        *result = sum;
    }
}

// residual = (idx == 0) ? 0 : rhs - ap. Cell 0 is forced to exactly 0
// regardless of rhs[0]'s value -- rhs[0] is never read here, matching the
// CPU original (StaggeredCavityBase3D::projectToDivergenceFree), where
// residual[0] is likewise forced to 0 independent of rhs[0].
__global__ void initResidualKernel(double* residual, const double* rhs, const double* ap, std::size_t n) {
    const std::size_t idx = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }
    residual[idx] = (idx == 0) ? 0.0 : rhs[idx] - ap[idx];
}

// p += alpha * direction, cell 0 untouched (stays exactly the caller's
// initial guess for that cell, for the whole solve).
__global__ void updatePressureKernel(double* p, const double* direction, double alpha, std::size_t n) {
    const std::size_t idx = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= n || idx == 0) {
        return;
    }
    p[idx] += alpha * direction[idx];
}

// residual -= alpha * ad, in place (idx == 0) forced to 0. Safe in place:
// each thread only ever reads and writes its own idx, and nothing later
// in the same iteration needs the pre-update residual value at this same
// index -- the CPU original's separate "newResidual" array is a stylistic
// choice there, not a true aliasing requirement.
__global__ void updateResidualKernel(double* residual, const double* ad, double alpha, std::size_t n) {
    const std::size_t idx = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }
    residual[idx] = (idx == 0) ? 0.0 : residual[idx] - alpha * ad[idx];
}

// direction = residual + beta * direction, in place (idx == 0) forced to
// 0. Also safe in place: `residual` here is *already* the freshly-updated
// value from updateResidualKernel (a separate, already-completed launch),
// and each thread reads the old direction[idx] once before overwriting
// it -- no cross-thread dependency.
__global__ void updateDirectionKernel(double* direction, const double* residual, double beta, std::size_t n) {
    const std::size_t idx = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (idx >= n) {
        return;
    }
    direction[idx] = (idx == 0) ? 0.0 : residual[idx] + beta * direction[idx];
}

} // namespace

struct ConjugateGradientSolverCuda::Impl {
    double* p = nullptr;
    double* rhs = nullptr;
    double* residual = nullptr;
    double* direction = nullptr;
    double* ad = nullptr;
    double* partialSums = nullptr;
    double* reducedScalar = nullptr;
    std::size_t numBlocks = 0;

    // Unconditional cudaFree on every pointer, whether it was ever
    // successfully allocated or not (still nullptr if construction failed
    // partway or no device was available) -- cudaFree(nullptr) is a
    // defined no-op, so this is safe in every case.
    ~Impl() {
        cudaFree(p);
        cudaFree(rhs);
        cudaFree(residual);
        cudaFree(direction);
        cudaFree(ad);
        cudaFree(partialSums);
        cudaFree(reducedScalar);
    }

    // Launches the two-pass reduction and returns the single reduced
    // scalar via one 8-byte device-to-host copy.
    double reduceDot(const double* a, const double* b, std::size_t n) const {
        dotProductPartialKernel<<<static_cast<unsigned>(numBlocks), kBlockSize1D>>>(a, b, partialSums, n);
        sumReduceKernel<<<1, kBlockSize1D>>>(partialSums, numBlocks, reducedScalar);
        double result = 0.0;
        cudaMemcpy(&result, reducedScalar, sizeof(double), cudaMemcpyDeviceToHost);
        return result;
    }
};

ConjugateGradientSolverCuda::ConjugateGradientSolverCuda(std::size_t nx, std::size_t ny, std::size_t nz,
                                                          double dx, double dy, double dz)
    : impl_(std::make_unique<Impl>()), nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz),
      available_(false) {
    int deviceCount = 0;
    if (cudaGetDeviceCount(&deviceCount) != cudaSuccess || deviceCount == 0) {
        return; // available_ stays false; impl_'s still-null buffers are freed harmlessly by ~Impl()
    }

    const std::size_t n = nx_ * ny_ * nz_;
    impl_->numBlocks = std::max<std::size_t>(
        1, std::min<std::size_t>((n + kBlockSize1D - 1) / kBlockSize1D, kMaxReductionBlocks));

    // Every allocation checked individually -- unlike PoissonOperatorCuda
    // (which allocates two transient buffers per call and never checks),
    // this class allocates seven persistent buffers once, so a large-grid
    // VRAM exhaustion is a real, not theoretical, failure mode. impl_ is
    // already a fully-constructed unique_ptr<Impl> at this point (built in
    // the member-initializer list above), so a throw here still unwinds
    // through ~Impl(), freeing whatever subset already succeeded.
    const auto checkedMalloc = [](double** ptr, std::size_t count) {
        const cudaError_t status = cudaMalloc(ptr, count * sizeof(double));
        if (status != cudaSuccess) {
            throw std::runtime_error(std::string("ConjugateGradientSolverCuda: cudaMalloc failed: ") +
                                      cudaGetErrorString(status));
        }
    };
    checkedMalloc(&impl_->p, n);
    checkedMalloc(&impl_->rhs, n);
    checkedMalloc(&impl_->residual, n);
    checkedMalloc(&impl_->direction, n);
    checkedMalloc(&impl_->ad, n);
    checkedMalloc(&impl_->partialSums, impl_->numBlocks);
    checkedMalloc(&impl_->reducedScalar, 1);

    available_ = true;
}

ConjugateGradientSolverCuda::~ConjugateGradientSolverCuda() = default;

ConjugateGradientSolverCuda::Result ConjugateGradientSolverCuda::solve(const std::vector<double>& rhs,
                                                                        const std::vector<double>& initialGuess,
                                                                        double tolerance,
                                                                        std::size_t maxIterations) const {
    if (!available_) {
        throw std::runtime_error("ConjugateGradientSolverCuda::solve: no CUDA-capable device available");
    }
    const std::size_t n = nx_ * ny_ * nz_;
    if (rhs.size() != n || initialGuess.size() != n) {
        throw std::invalid_argument("ConjugateGradientSolverCuda::solve: rhs/initialGuess size mismatch");
    }
    if (maxIterations == 0) {
        maxIterations = n;
    }

    Impl& impl = *impl_;
    cudaMemcpy(impl.p, initialGuess.data(), n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(impl.rhs, rhs.data(), n * sizeof(double), cudaMemcpyHostToDevice);

    const double ax = 1.0 / (dx_ * dx_);
    const double ay = 1.0 / (dy_ * dy_);
    const double az = 1.0 / (dz_ * dz_);
    const double weightTotal = 2.0 * ax + 2.0 * ay + 2.0 * az;

    const dim3 blockDim3(8, 8, 8);
    const dim3 gridDim3(static_cast<unsigned>((nx_ + 7) / 8), static_cast<unsigned>((ny_ + 7) / 8),
                         static_cast<unsigned>((nz_ + 7) / 8));
    const unsigned gridDim1D = static_cast<unsigned>((n + kBlockSize1D - 1) / kBlockSize1D);

    // Ap = applyLaplacian(p), once, before the loop -- mirrors
    // projectToDivergenceFree()'s own pre-loop residual initialization.
    applyLaplacianKernel<<<gridDim3, blockDim3>>>(impl.p, impl.ad, nx_, ny_, nz_, ax, ay, az, weightTotal);
    initResidualKernel<<<gridDim1D, kBlockSize1D>>>(impl.residual, impl.rhs, impl.ad, n);
    cudaMemcpy(impl.direction, impl.residual, n * sizeof(double), cudaMemcpyDeviceToDevice);

    double residualDotResidual = impl.reduceDot(impl.residual, impl.residual, n);

    std::size_t iteration = 0;
    bool brokeDown = false;
    for (; iteration < maxIterations; ++iteration) {
        if (std::sqrt(residualDotResidual) < tolerance) {
            break;
        }
        applyLaplacianKernel<<<gridDim3, blockDim3>>>(impl.direction, impl.ad, nx_, ny_, nz_, ax, ay, az,
                                                       weightTotal);
        const double directionDotAd = impl.reduceDot(impl.direction, impl.ad, n);
        if (directionDotAd == 0.0) {
            brokeDown = true;
            break;
        }
        const double alpha = residualDotResidual / directionDotAd;
        updatePressureKernel<<<gridDim1D, kBlockSize1D>>>(impl.p, impl.direction, alpha, n);
        updateResidualKernel<<<gridDim1D, kBlockSize1D>>>(impl.residual, impl.ad, alpha, n);
        const double newResidualDotResidual = impl.reduceDot(impl.residual, impl.residual, n);
        const double beta = newResidualDotResidual / residualDotResidual;
        updateDirectionKernel<<<gridDim1D, kBlockSize1D>>>(impl.direction, impl.residual, beta, n);
        residualDotResidual = newResidualDotResidual;
    }

    Result result;
    result.pressure.resize(n);
    cudaDeviceSynchronize();
    cudaMemcpy(result.pressure.data(), impl.p, n * sizeof(double), cudaMemcpyDeviceToHost);

    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("ConjugateGradientSolverCuda::solve: CUDA error: ") +
                                  cudaGetErrorString(status));
    }

    result.iterations = iteration;
    result.residualNorm = std::sqrt(residualDotResidual);
    result.brokeDown = brokeDown;
    result.converged = result.residualNorm < tolerance;
    return result;
}

} // namespace aether::gpu
