#pragma once

#include <cstddef>
#include <memory>
#include <vector>

namespace aether::gpu {

// Fase 4.1 of ROADMAP.md: the same matrix-free conjugate-gradient
// pressure-projection solve every StaggeredCavityBase3D-derived 3D solver
// runs on the CPU (see StaggeredCavityBase3D::projectToDivergenceFree),
// but with the entire iterative loop -- including both dot-product
// reductions per iteration -- resident on the GPU. PoissonOperatorCuda
// copies host<->device on every apply() call, which is fine for
// validating the stencil kernel alone but means a CG solve (which calls
// that operator dozens of times per pressure projection) never actually
// gets faster: the PCIe round trip dominates every iteration. This class
// transfers the initial rhs/guess once and the final pressure once; every
// intermediate field (p, residual, direction, the operator's own output)
// stays on the device for the whole solve. Only two 8-byte scalars (the
// two reduced dot products) cross the bus per iteration, to compute
// alpha/beta and check convergence on the host.
//
// **Scope, deliberate**: this class is standalone and untested against
// engine/solver -- it is not yet wired into StaggeredCavityBase3D's own
// projectToDivergenceFree(). That integration is a separate, later step
// (same reasoning DelaunayTetrahedralization3D's O(N^2) fix used to leave
// insertSteinerPoint() untouched: a narrower, independently-verifiable
// increment is faster to review and verify completely).
//
// **Validation is not bit-for-bit against the CPU**, unlike
// PoissonOperatorCuda's own stencil-only comparison -- and deliberately
// so: a GPU dot-product reduction sums in a different (but fixed,
// reproducible run-to-run -- see the .cu file's own comment on why
// atomics were rejected) order than the CPU's sequential loop, so the
// two trajectories are expected to diverge slightly after enough
// iterations. What's checked instead (see gpu_tests.cpp): the GPU's own
// returned pressure field independently re-satisfies residual = rhs - A*p
// on the CPU within a measured tolerance, and the pinned cell (index 0)
// matches the initial guess bit-for-bit unconditionally (a cheap,
// always-applicable regression check for the most likely bug: one of the
// three update kernels forgetting its own idx==0 guard).
class ConjugateGradientSolverCuda {
public:
    // nx/ny/nz: grid dimensions (cell index i + j*nx + k*nx*ny, same
    // layout every staggered 3D solver and PoissonOperatorCuda use).
    // dx/dy/dz: cell spacing along each axis. Device buffers sized to
    // nx*ny*nz are allocated here, eagerly, and reused by every solve()
    // call on this instance -- never reallocated per call.
    ConjugateGradientSolverCuda(std::size_t nx, std::size_t ny, std::size_t nz, double dx, double dy,
                                 double dz);
    ~ConjugateGradientSolverCuda();

    ConjugateGradientSolverCuda(const ConjugateGradientSolverCuda&) = delete;
    ConjugateGradientSolverCuda& operator=(const ConjugateGradientSolverCuda&) = delete;

    struct Result {
        std::vector<double> pressure; // size nx*ny*nz
        std::size_t iterations = 0;
        double residualNorm = 0.0; // sqrt(residualDotResidual) at loop exit, as the GPU itself computed it
        bool converged = false;    // residualNorm < tolerance
        bool brokeDown = false;    // direction . (A * direction) == 0.0 exactly (CG breakdown)
    };

    // Mirrors StaggeredCavityBase3D::projectToDivergenceFree()'s own CG
    // loop exactly: cell 0 is pinned (removes the pure-Neumann system's
    // null space) -- rhs[0]'s value is never read, since residual[0] is
    // forced to 0 regardless, the same as the CPU version. maxIterations
    // == 0 means "use nx*ny*nz", the same cap the CPU loop uses (can't be
    // a default argument directly: n isn't known at header-parse time).
    // Throws std::runtime_error if no CUDA-capable device is available
    // (check available() first to avoid the exception), and
    // std::invalid_argument if rhs/initialGuess don't both have size
    // nx*ny*nz.
    Result solve(const std::vector<double>& rhs, const std::vector<double>& initialGuess,
                 double tolerance = 1e-10, std::size_t maxIterations = 0) const;

    // True if a CUDA-capable device was found and every device buffer
    // this instance needs was successfully allocated at construction --
    // check this before calling solve() if the caller wants to fall back
    // to a CPU path instead of catching an exception.
    bool available() const { return available_; }

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    std::size_t nx_;
    std::size_t ny_;
    std::size_t nz_;
    double dx_;
    double dy_;
    double dz_;
    bool available_;
};

} // namespace aether::gpu
