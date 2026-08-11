#pragma once

#include <cstddef>
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

protected:
    StaggeredCavityBase3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX, double lengthY,
                           double lengthZ, double viscosity, double lidVelocity);
    ~StaggeredCavityBase3D() = default;

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

private:
    const std::vector<double>* eddyViscosity_ = nullptr;
};

} // namespace aether::solver
