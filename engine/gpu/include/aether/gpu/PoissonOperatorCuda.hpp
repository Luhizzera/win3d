#pragma once

#include <cstddef>
#include <vector>

namespace aether::gpu {

// Module 10, first pass: the project's first CUDA kernel.
//
// **What it does**: applies the exact same matrix-free discrete Poisson
// operator every pressure-projection solver in this codebase already runs
// on the CPU, dozens of times per call to solveConjugateGradient() /
// projectToDivergenceFree() -- see StaggeredCavityBase3D::applyLaplacian
// (and its several standalone 2D/3D siblings predating that class) for the
// CPU original this mirrors line-for-line: the standard 7-point
// finite-volume Laplacian, "clamp to the boundary cell" Neumann treatment
// on all six faces, and cell 0 pinned to remove the pure-Neumann system's
// null space.
//
// **Why this operator, first**: it is the single most reused piece of
// numerics in the codebase (every explicit Navier-Stokes solver here, 2D
// or 3D, laminar or any of the five turbulence closures, calls something
// of exactly this shape once per pressure-projection step), and it is
// embarrassingly parallel -- every output cell depends on exactly 6 fixed
// neighbours with no cross-cell accumulation order to get subtly wrong,
// which makes it as close to a risk-free first GPU port as this engine
// has.
//
// **The validation this enables, and why it's stronger than usual**: MSVC
// and nvcc are both IEEE-754-strict by default (neither enables fast-math
// unless explicitly asked, and this class asks for neither -std=fast-math
// nor --use_fast_math), and this kernel performs the identical sequence of
// additions/multiplications the CPU version does, in the identical order,
// on identical binary64 doubles. So GPU and CPU outputs are compared for
// **exact bit-for-bit agreement**, not "close enough" -- a strictly
// stronger check than this project's usual tolerance-based validation,
// possible here specifically because a plain stencil (unlike a reduction
// such as a dot product) has no summation-order ambiguity to begin with.
//
// **Deliberately not yet**: a persistent device-resident solve (this
// class copies host->device, launches, copies back -- once per apply()
// call, no state kept on the device between calls). Fine for validating
// the kernel itself; a real speedup needs the whole CG loop (dot products
// included) resident on the GPU so the PCIe transfer isn't paid every
// iteration. That is explicitly future Module 10 work, not attempted here
// so the first CUDA code in this project stays small enough to review and
// validate completely.
class PoissonOperatorCuda {
public:
    // nx/ny/nz: grid dimensions (cell index i + j*nx + k*nx*ny, the same
    // layout StructuredGrid3D and every staggered 3D solver here uses).
    // dx/dy/dz: cell spacing along each axis.
    PoissonOperatorCuda(std::size_t nx, std::size_t ny, std::size_t nz, double dx, double dy, double dz);
    ~PoissonOperatorCuda();

    PoissonOperatorCuda(const PoissonOperatorCuda&) = delete;
    PoissonOperatorCuda& operator=(const PoissonOperatorCuda&) = delete;

    // Applies the operator to x (size nx*ny*nz). Throws std::runtime_error
    // if no CUDA-capable device was found at construction, and
    // std::invalid_argument if x.size() doesn't match nx*ny*nz.
    std::vector<double> apply(const std::vector<double>& x) const;

    // True if a CUDA-capable device was found at construction time --
    // check this before calling apply() if the caller wants to fall back
    // to a CPU path instead of catching an exception.
    bool available() const { return available_; }

private:
    std::size_t nx_;
    std::size_t ny_;
    std::size_t nz_;
    double dx_;
    double dy_;
    double dz_;
    bool available_;
};

} // namespace aether::gpu
