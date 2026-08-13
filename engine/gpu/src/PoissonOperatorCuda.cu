#include "aether/gpu/PoissonOperatorCuda.hpp"

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>

namespace aether::gpu {

namespace {

// Line-for-line the same stencil as StaggeredCavityBase3D::applyLaplacian
// on the CPU: for the pinned cell (idx==0) the operator is the identity;
// otherwise each of the six neighbour lookups clamps to the current cell's
// own layer at a domain boundary (the "ghost = interior" Neumann mirror),
// and the result is weightTotal*x[idx] minus the weighted neighbour sum.
// One thread per cell -- nx*ny*nz is small enough at every resolution this
// project uses (up to a few hundred thousand cells) that one cell per
// thread, no shared-memory tiling, is the right amount of complexity for a
// first kernel whose entire purpose is to be checked bit-for-bit against
// the CPU original.
__global__ void poissonKernel(const double* x, double* result, std::size_t nx, std::size_t ny, std::size_t nz,
                               double ax, double ay, double az, double weightTotal) {
    const auto i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto j = static_cast<std::size_t>(blockIdx.y) * blockDim.y + threadIdx.y;
    const auto k = static_cast<std::size_t>(blockIdx.z) * blockDim.z + threadIdx.z;
    if (i >= nx || j >= ny || k >= nz) {
        return;
    }
    const std::size_t idx = i + j * nx + k * nx * ny;
    if (idx == 0) {
        result[idx] = x[idx];
        return;
    }
    const std::size_t left = (i > 0 ? i - 1 : 0) + j * nx + k * nx * ny;
    const std::size_t right = (i + 1 < nx ? i + 1 : nx - 1) + j * nx + k * nx * ny;
    const std::size_t down = i + (j > 0 ? j - 1 : 0) * nx + k * nx * ny;
    const std::size_t up = i + (j + 1 < ny ? j + 1 : ny - 1) * nx + k * nx * ny;
    const std::size_t back = i + j * nx + (k > 0 ? k - 1 : 0) * nx * ny;
    const std::size_t front = i + j * nx + (k + 1 < nz ? k + 1 : nz - 1) * nx * ny;

    const double weightedSum = ax * (x[left] + x[right]) + ay * (x[down] + x[up]) + az * (x[back] + x[front]);
    result[idx] = weightTotal * x[idx] - weightedSum;
}

} // namespace

PoissonOperatorCuda::PoissonOperatorCuda(std::size_t nx, std::size_t ny, std::size_t nz, double dx, double dy,
                                          double dz)
    : nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz), available_(false) {
    int deviceCount = 0;
    available_ = cudaGetDeviceCount(&deviceCount) == cudaSuccess && deviceCount > 0;
}

PoissonOperatorCuda::~PoissonOperatorCuda() = default;

std::vector<double> PoissonOperatorCuda::apply(const std::vector<double>& x) const {
    if (!available_) {
        throw std::runtime_error("PoissonOperatorCuda::apply: no CUDA-capable device available");
    }
    const std::size_t n = nx_ * ny_ * nz_;
    if (x.size() != n) {
        throw std::invalid_argument("PoissonOperatorCuda::apply: x.size() does not match nx*ny*nz");
    }

    double* deviceX = nullptr;
    double* deviceResult = nullptr;
    cudaMalloc(&deviceX, n * sizeof(double));
    cudaMalloc(&deviceResult, n * sizeof(double));
    cudaMemcpy(deviceX, x.data(), n * sizeof(double), cudaMemcpyHostToDevice);

    const double ax = 1.0 / (dx_ * dx_);
    const double ay = 1.0 / (dy_ * dy_);
    const double az = 1.0 / (dz_ * dz_);
    const double weightTotal = 2.0 * ax + 2.0 * ay + 2.0 * az;

    const dim3 blockDim(8, 8, 8);
    const dim3 gridDim(static_cast<unsigned>((nx_ + 7) / 8), static_cast<unsigned>((ny_ + 7) / 8),
                        static_cast<unsigned>((nz_ + 7) / 8));
    poissonKernel<<<gridDim, blockDim>>>(deviceX, deviceResult, nx_, ny_, nz_, ax, ay, az, weightTotal);
    cudaDeviceSynchronize();

    std::vector<double> result(n);
    cudaMemcpy(result.data(), deviceResult, n * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(deviceX);
    cudaFree(deviceResult);

    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("PoissonOperatorCuda::apply: CUDA error: ") +
                                  cudaGetErrorString(status));
    }
    return result;
}

} // namespace aether::gpu
