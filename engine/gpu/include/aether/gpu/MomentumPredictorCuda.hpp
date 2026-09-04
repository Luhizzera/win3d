#pragma once

#include <cstddef>
#include <memory>
#include <vector>

namespace aether::gpu {

// Fase 4.2 of ROADMAP.md: StaggeredCavityBase3D::computeMomentumPredictor
// (engine/solver/src/StaggeredCavityBase3D.cpp), ported to CUDA -- the
// explicit predictor step every one of the six StaggeredCavityBase3D-
// derived 3D solvers runs once per time step, before the pressure
// projection Fase 4.1's ConjugateGradientSolverCuda already covers.
//
// **Scope, deliberate, verified against the source rather than assumed**:
// only ConvectionScheme::Central is ported -- the only scheme any of the
// six derived closures actually constructs with today. For Central,
// StaggeredCavityBase3D::schemeTransportValue's own body returns its
// `centralValue` argument completely unchanged, never touching the
// upwind/downwind/limiter machinery at all -- so every far-point lookup
// that exists solely to feed that limiter is dead computation for the
// case this class covers, and dropping it is behavior-preserving, not an
// approximation. Turbulent viscosity (nu_t) IS fully supported, via an
// optional per-cell `nut` field -- unlike the convection scheme, every
// one of the five turbulent closures needs this to produce correct
// output.
//
// **Standalone, not yet wired into engine/solver** -- same reasoning as
// ConjugateGradientSolverCuda: a narrower, independently-verifiable
// increment is faster to review and verify completely. Real integration
// (and the real speedup measurement that only means something once this
// runs inside an actual time-stepping loop) is a separate, later step.
//
// **Validation is bit-for-bit against the CPU**, unlike
// ConjugateGradientSolverCuda's own tolerance-based comparison -- and
// correctly so: this kernel has no reduction anywhere (every output cell
// is a fixed, order-independent arithmetic expression over neighbor
// reads), so with -fmad=false already applying project-wide (see
// engine/gpu/CMakeLists.txt), GPU and CPU are expected to agree exactly,
// the same bar PoissonOperatorCuda's own test already uses.
class MomentumPredictorCuda {
public:
    // Mirrors StaggeredCavityBase3D's own protected constructor
    // parameters exactly.
    MomentumPredictorCuda(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX, double lengthY,
                           double lengthZ, double viscosity, double lidVelocity);
    ~MomentumPredictorCuda();

    MomentumPredictorCuda(const MomentumPredictorCuda&) = delete;
    MomentumPredictorCuda& operator=(const MomentumPredictorCuda&) = delete;

    struct Result {
        std::vector<double> uStar; // size (nx+1)*ny*nz
        std::vector<double> vStar; // size nx*(ny+1)*nz
        std::vector<double> wStar; // size nx*ny*(nz+1)
    };

    // u/v/w must match StaggeredCavityBase3D's own staggered sizes
    // exactly: (nx+1)*ny*nz, nx*(ny+1)*nz, nx*ny*(nz+1) respectively.
    // `nut` is either empty (laminar -- the device nut buffer is
    // explicitly re-zeroed for this call, since it is reused across
    // calls and a prior turbulent call could otherwise leave stale
    // values behind) or exactly nx*ny*nz, cell-centered, same convention
    // as StaggeredCavityBase3D::nutAt's own registered field. Boundary
    // faces of uStar/vStar/wStar (i in {0,nx} for u, j in {0,ny} for v,
    // k in {0,nz} for w) are copied from u/v/w verbatim, never computed
    // -- matching the CPU contract that computeMomentumPredictor only
    // ever writes interior faces. Throws std::runtime_error if no
    // CUDA-capable device is available, std::invalid_argument on a size
    // mismatch.
    Result predict(const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& w,
                    const std::vector<double>& nut, double dt) const;

    // True if a CUDA-capable device was found and every device buffer
    // this instance needs was successfully allocated at construction.
    bool available() const { return available_; }

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    std::size_t nx_;
    std::size_t ny_;
    std::size_t nz_;
    double lengthX_;
    double lengthY_;
    double lengthZ_;
    double dx_;
    double dy_;
    double dz_;
    double viscosity_;
    double lidVelocity_;
    bool available_;
};

} // namespace aether::gpu
