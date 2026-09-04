#pragma once

#include <cstddef>
#include <memory>
#include <vector>

namespace aether::solver {

// Shared machinery for every 3D lid-driven cavity solver in this project,
// extracted in the cleanup pass before the roadmap's product modules
// (9-14) began.
//
// **Why it was extracted, and why then**: six classes
// (StaggeredLidDrivenCavitySolver3D plus the mixing-length, k-epsilon,
// k-omega SST, Smagorinsky LES and SST-DES closures) had each grown their
// own verbatim copy of the staggered/MAC grid plumbing -- roughly 300
// lines apiece of indexing, ghost-cell mirrors, the pressure projection,
// and the momentum predictor. That was tolerable while closures were being
// added one at a time, but modules 10 (GPU) and 9/11/13 (UI, persistence,
// API) all bind against exactly this surface: porting the *same* kernel
// six times, or six-way-updating a shared API, is the expensive version of
// this problem. Cheaper to pay it here.
//
// **The extraction was verified to be behavior-preserving before it was
// made, not assumed.** Every function below was compared across all six
// classes by normalizing whitespace/comments and hashing the bodies:
// uAt, vAt, wAt, maxDivergence, applyLaplacian, projectToDivergenceFree,
// lidVelocityAt, wallDistanceAt and dot were **byte-identical** in all six;
// the momentum predictor was byte-identical across all five turbulent
// closures. Only two things differed, both understood:
//   - pAt: the LES class spelled the clamp with std::clamp and the others
//     with explicit if/else -- semantically the same; the std::clamp form
//     is kept here.
//   - stableTimeStep: the laminar class has no nu_t to add. Handled by
//     maxEddyViscosity() returning exactly 0.0 when no eddy-viscosity
//     field is registered, which makes `viscosity_ + 0.0` reproduce the
//     original expression bit-for-bit.
//
// **No virtual functions, deliberately.** The obvious design -- a virtual
// eddyViscosityAt() hook -- would put a virtual dispatch in the innermost
// momentum loop, called six times per cell per velocity component. Instead
// a derived closure registers a pointer to its own nu_t vector via
// setEddyViscosityField(); nutAt() is a plain inlinable function that
// returns 0.0 when that pointer is null (the laminar case). One predictable
// branch instead of a vtable lookup, and the class stays non-polymorphic
// (hence the protected non-virtual destructor: deleting through a base
// pointer is neither needed nor allowed here).
class StaggeredCavityBase3D {
public:
    // How the convected velocity is reconstructed at a face, mirroring
    // LidDrivenCavitySolver2D::ConvectionScheme (DIVIDA_TECNICA.md 4.4) --
    // not the same type, because sharing it would mean the 2D class's
    // enum stops being self-contained for no benefit here, but the same
    // three points on the same line:
    //
    //   phi = phi_upwind + psi * (phi_central - phi_upwind)
    //
    // Central is the default because it is what every one of the six
    // classes built on this base has always used: every published number
    // stays reproducible bit-for-bit until a caller opts into one of the
    // other two. See computeMomentumPredictor() for why "central" here
    // means something slightly different from the 2D collocated case: this
    // grid is staggered, so the value that already plays "the central
    // estimate" in the original conservative formula (a plain average
    // between two neighboring staggered points) is reused unchanged as
    // Central's output, and only the *other* two schemes reconstruct a
    // different value at that same location.
    enum class ConvectionScheme {
        Central,
        FirstOrderUpwind,
        LimitedLinearUpwind,
    };

    double u(std::size_t i, std::size_t j, std::size_t k) const { return u_[indexU(i, j, k)]; }
    double v(std::size_t i, std::size_t j, std::size_t k) const { return v_[indexV(i, j, k)]; }
    double w(std::size_t i, std::size_t j, std::size_t k) const { return w_[indexW(i, j, k)]; }
    double pressure(std::size_t i, std::size_t j, std::size_t k) const { return p_[indexP(i, j, k)]; }
    double time() const { return time_; }

    // Conservative explicit-stepping limit: the diffusive bound uses
    // molecular plus (if a closure registered one) the current maximum
    // eddy viscosity, the convective bound uses the lid speed.
    double stableTimeStep() const;

    // Max |div(u)| over all cells -- the mass-conservation diagnostic every
    // one of these solvers' tests checks every step.
    double maxDivergence() const;

    // True iff the GPU path (Fase 4 of ROADMAP.md) is actually live for
    // this instance: `useGpu=true` was requested at construction AND a
    // CUDA-capable device was found AND the device buffers it needs were
    // successfully allocated. False whenever any of those isn't true --
    // including a build compiled without CUDA support at all -- so a
    // caller never needs to know about AETHER_HAVE_GPU to use this
    // correctly; it behaves identically regardless of build configuration.
    // Public (not protected) for the same reason u()/v()/w()/... are:
    // automatically inherited as public by every derived class with no
    // `using` declaration needed, unlike loadState().
    bool gpuActive() const { return gpuActive_; }

protected:
    // Overwrites u/v/w/p and the simulated time with externally supplied
    // data -- the counterpart to reading them out to save a checkpoint (see
    // engine/persistence/FieldArchive). Protected, not public: only the
    // *laminar* closure (no extra turbulence fields) can safely expose this
    // as-is via a `using` declaration -- the five turbulent closures each
    // carry additional state (k, epsilon, omega, nut, ...) this method
    // knows nothing about, so exposing it on them directly would silently
    // resume only part of their state. Throws std::invalid_argument if any
    // field's size doesn't match the grid's own u_/v_/w_/p_ sizes.
    void loadState(std::vector<double> u, std::vector<double> v, std::vector<double> w, std::vector<double> p,
                   double time);

    // useGpu: opt-in (Fase 4 of ROADMAP.md), defaults to false -- every
    // existing derived-class constructor passes its base-constructor
    // arguments positionally, never by name, so this new trailing
    // parameter cannot change behavior for any of them. When true but the
    // GPU path can't actually be brought up (no CUDA support compiled in,
    // no device found at runtime, or device allocation failed), the
    // constructor falls back to the CPU path silently except for a
    // diagnostic on stderr -- see the .cpp for why that's a fallback
    // rather than a thrown exception (an unavailable GPU is a runtime
    // resource condition, not a caller mistake, unlike this class's other
    // validation throws).
    StaggeredCavityBase3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX, double lengthY,
                           double lengthZ, double viscosity, double lidVelocity,
                           ConvectionScheme convection = ConvectionScheme::Central, bool useGpu = false);
    // Declaration only -- GpuState (below) is incomplete here. Defined in
    // the .cpp, the one file that ever needs it complete.
    ~StaggeredCavityBase3D();

    std::size_t indexU(std::size_t i, std::size_t j, std::size_t k) const {
        return i + j * (nx_ + 1) + k * (nx_ + 1) * ny_;
    }
    std::size_t indexV(std::size_t i, std::size_t j, std::size_t k) const {
        return i + j * nx_ + k * nx_ * (ny_ + 1);
    }
    std::size_t indexW(std::size_t i, std::size_t j, std::size_t k) const { return i + j * nx_ + k * nx_ * ny_; }
    std::size_t indexP(std::size_t i, std::size_t j, std::size_t k) const { return i + j * nx_ + k * nx_ * ny_; }

    // Lid tapered as sin^2 in both directions so it vanishes at the walls
    // it meets -- the regularization that removes the corner singularity a
    // discontinuous lid would otherwise impose (see LidDrivenCavitySolver2D
    // for where that was first diagnosed).
    double lidVelocityAt(double x, double y) const;
    // Distance from a cell center to the nearest of the six walls.
    double wallDistanceAt(std::size_t i, std::size_t j, std::size_t k) const;

    // Ghost-mirrored velocity/pressure lookups: no-slip on five walls, the
    // tapered lid on z = lengthZ, zero-gradient for pressure everywhere.
    double uAt(long long i, long long j, long long k) const;
    double vAt(long long i, long long j, long long k) const;
    double wAt(long long i, long long j, long long k) const;
    double pAt(long long i, long long j, long long k) const;

    // Registers a closure's own nu_t field (cell-centered, indexP layout).
    // Not owned; the vector must outlive this object and must not be
    // reallocated afterwards. Leaving it unregistered means laminar flow.
    void setEddyViscosityField(const std::vector<double>* field) { eddyViscosity_ = field; }

    // nu_t at a cell, 0.0 outside the domain or when no field is
    // registered -- matching what each closure's own nutAt() did before
    // this base existed.
    double nutAt(long long i, long long j, long long k) const;
    double maxEddyViscosity() const;

    // Explicit predictor for all three momentum components: conservative
    // convection plus face-weighted (nu + nu_t) diffusion, on the interior
    // faces only. Writes uStar/vStar/wStar, which must already be sized
    // like u_/v_/w_ (callers pass copies of them, so boundary faces keep
    // their prescribed values).
    void computeMomentumPredictor(std::vector<double>& uStar, std::vector<double>& vStar,
                                   std::vector<double>& wStar, double dt) const;

    // Reconstructs the value transported through a location shared by two
    // neighboring staggered points, given the plain average of those two
    // points (`centralValue`, always computed the same way regardless of
    // scheme) and the velocity whose sign decides which of the two is
    // upwind (`convectingVelocity` -- the same quantity as centralValue for
    // a component's self-convection term, a different one for its cross
    // terms). `near0`/`near1` are the two points centralValue averages;
    // `far0`/`far1` are one step beyond each, in the same order, for the
    // limiter's ratio -- the caller is responsible for clamping those to a
    // valid index when the true neighbor would fall outside the domain
    // (see computeMomentumPredictor.cpp for why that only matters for a
    // component's own staggered direction).
    //
    // Central returns centralValue unchanged -- not merely equal to it,
    // literally the same double, no arithmetic performed -- which is what
    // makes ConvectionScheme::Central bit-for-bit identical to this
    // predictor's original (pre-4.4) formula.
    double schemeTransportValue(double centralValue, double convectingVelocity, double near0, double near1,
                                 double far0, double far1) const;

    // Matrix-free Poisson operator for the pressure correction, with cell 0
    // pinned to remove the pure-Neumann null space.
    std::vector<double> applyLaplacian(const std::vector<double>& x) const;
    static double dot(const std::vector<double>& a, const std::vector<double>& b);

    // Chorin projection: solve for the pressure that makes the predicted
    // velocity divergence-free (matrix-free CG), then correct u/v/w.
    void projectToDivergenceFree(std::vector<double>& uStar, std::vector<double>& vStar,
                                  std::vector<double>& wStar, double dt);

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
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> w_;
    std::vector<double> p_;
    double time_ = 0.0;
    ConvectionScheme convection_;

private:
    const std::vector<double>* eddyViscosity_ = nullptr;

    // Forward-declared, incomplete here -- the same pImpl idiom
    // aether::gpu::ConjugateGradientSolverCuda/MomentumPredictorCuda's own
    // headers already use for their own `Impl`. The header's declared
    // layout is therefore a single pointer, byte-identical in every
    // translation unit regardless of whether AETHER_HAVE_GPU is ever
    // defined anywhere -- no translation unit can ever disagree about
    // sizeof(StaggeredCavityBase3D). GpuState's real definition (or an
    // empty stand-in when compiled without CUDA support) lives entirely in
    // StaggeredCavityBase3D.cpp, the only file that needs to know which.
    struct GpuState;
    std::unique_ptr<GpuState> gpuState_;
    bool gpuActive_ = false;
};

} // namespace aether::solver
