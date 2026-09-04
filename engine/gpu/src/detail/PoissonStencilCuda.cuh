#pragma once

#include <cstddef>

namespace aether::gpu::detail {

// The one piece of index arithmetic PoissonOperatorCuda.cu and
// ConjugateGradientSolverCuda.cu must never be allowed to drift apart on:
// the 7-point finite-volume Laplacian stencil, cell 0 pinned to identity
// (removes the pure-Neumann system's null space), the other six neighbour
// lookups clamped to the current cell's own layer at a domain boundary
// (the "ghost = interior" Neumann mirror). Line-for-line the same formula
// as StaggeredCavityBase3D::applyLaplacian on the CPU (engine/solver/src/
// StaggeredCavityBase3D.cpp) and, before this extraction, as
// PoissonOperatorCuda.cu's own poissonKernel had inlined directly.
//
// A __device__ function in a header, not a __global__ kernel: each caller
// keeps its own launch configuration and grid-stride/1D indexing scheme,
// only the per-cell formula itself is shared. __forceinline__ so this
// costs nothing at the call site relative to the old fully-inlined
// version -- verified by re-running PoissonOperatorCuda's own existing
// bit-exact test unchanged after this extraction.
__device__ __forceinline__ double poissonStencilValue(const double* x, std::size_t i, std::size_t j,
                                                        std::size_t k, std::size_t nx, std::size_t ny,
                                                        std::size_t nz, double ax, double ay, double az,
                                                        double weightTotal) {
    const std::size_t idx = i + j * nx + k * nx * ny;
    if (idx == 0) {
        return x[idx];
    }
    const std::size_t left = (i > 0 ? i - 1 : 0) + j * nx + k * nx * ny;
    const std::size_t right = (i + 1 < nx ? i + 1 : nx - 1) + j * nx + k * nx * ny;
    const std::size_t down = i + (j > 0 ? j - 1 : 0) * nx + k * nx * ny;
    const std::size_t up = i + (j + 1 < ny ? j + 1 : ny - 1) * nx + k * nx * ny;
    const std::size_t back = i + j * nx + (k > 0 ? k - 1 : 0) * nx * ny;
    const std::size_t front = i + j * nx + (k + 1 < nz ? k + 1 : nz - 1) * nx * ny;

    const double weightedSum = ax * (x[left] + x[right]) + ay * (x[down] + x[up]) + az * (x[back] + x[front]);
    return weightTotal * x[idx] - weightedSum;
}

} // namespace aether::gpu::detail
