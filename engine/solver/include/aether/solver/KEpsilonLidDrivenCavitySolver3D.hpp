#pragma once

#include "aether/solver/StaggeredCavityBase3D.hpp"

#include <cstddef>
#include <vector>

namespace aether::solver {

// The two-equation k-epsilon closure, extended from the real 2D cavity
// (KEpsilonLidDrivenCavitySolver2D) to the real 3D solid-wall cavity
// (StaggeredLidDrivenCavitySolver3D's MAC grid), mirroring the same step
// MixingLengthLidDrivenCavitySolver3D already took for the simpler
// algebraic closure. k, epsilon, nu_t and the production term all live at
// cell centers -- the same points pressure occupies -- so unlike the
// momentum equations (genuinely staggered, needing the own-axis/transverse
// treatment MixingLengthLidDrivenCavitySolver3D documents), k/epsilon
// transport is discretized on a plain colocated 6-neighbor stencil (two
// neighbors per axis, three axes), the direct 3D generalization of the 2D
// class's own 4-neighbor colocated treatment.
//
// **Wall treatment, momentum diffusion, warm-start and Patankar
// linearization are all carried over unchanged in spirit from
// KEpsilonLidDrivenCavitySolver2D** (see that class's own header comment
// for the full reasoning): epsilon's wall ghost uses the self-derived
// low-Re asymptotic relation epsilon_wall = 2*nu*k_nearWallCell/y_half^2
// (now evaluated against whichever of the three axes the crossed wall is
// on), k is fixed to exactly 0 at all six walls, velocity is warm-started
// from MixingLengthLidDrivenCavitySolver3D (same cold-start-deadlock
// avoidance), and epsilon's quadratic destruction term uses the same
// frozen-factor implicit form. The momentum equations reuse
// MixingLengthLidDrivenCavitySolver3D's own diffusion treatment verbatim
// (exact face-varying effective viscosity along a component's own axis,
// locally-averaged along the two transverse axes).
// The staggered-grid plumbing (indexing, ghost mirrors, momentum
// predictor, pressure projection, stableTimeStep, maxDivergence) lives in
// StaggeredCavityBase3D; only the closure itself is here.
class KEpsilonLidDrivenCavitySolver3D : public StaggeredCavityBase3D {
public:
    KEpsilonLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                     double lengthY, double lengthZ, double viscosity, double lidVelocity);

    void step(double dt);

    double k(std::size_t i, std::size_t j, std::size_t k) const;
    double epsilon(std::size_t i, std::size_t j, std::size_t k) const;
    double eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const;

private:
    // Cell-centered k; homogeneous Dirichlet ghost mirror (k=0 at any wall).
    double kAt(long long i, long long j, long long k) const;
    // Cell-centered epsilon; ghost mirrors the computed wall value
    // 2*nu*k_nearWallCell/y_half^2 -- see the class comment.
    double epsilonGhostAt(std::size_t i, std::size_t j, std::size_t k, int di, int dj, int dk) const;

    void updateEddyViscosityAndProduction();

    std::vector<double> k_;
    std::vector<double> epsilon_;
    std::vector<double> nut_;
    std::vector<double> production_;
};

} // namespace aether::solver
