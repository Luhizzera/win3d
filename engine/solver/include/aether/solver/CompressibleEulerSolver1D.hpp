#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Compressible Euler equations in 1D -- ROADMAP Fase 6.2, and a genuine
// regime change rather than an extension of anything upstream: every other
// solver in this project is incompressible (density fixed, pressure solved
// from a projection), and shocks need a hyperbolic system with its own
// numerical flux instead. This is the first slice, deliberately the
// simplest pair that still captures a real shock without oscillating --
// the same reasoning that took this project through first-order upwind
// before limited-linear-upwind, and mixing length before k-epsilon.
// Viscosity is left out on purpose: it already exists everywhere else in
// this project, and what is missing here is compressibility and shock
// capture, which Euler isolates cleanly.
//
// **Governing equations**, conservative form: dU/dt + dF(U)/dx = 0, with
// U = [rho, rho*u, E] and F(U) = [rho*u, rho*u^2+p, u*(E+p)], closed by the
// ideal-gas equation of state p = (gamma-1)*(E - 0.5*rho*u^2).
//
// **Numerical flux: Rusanov (local Lax-Friedrichs).**
//   F_face = 0.5*(F(U_L)+F(U_R)) - 0.5*s_max*(U_R-U_L)
// with s_max = max(|u_L|+c_L, |u_R|+c_R) and c = sqrt(gamma*p/rho) the
// local sound speed. The simplest robust flux this system has: no
// eigenvalue decomposition (Roe) and no wave-speed splitting (HLLC), just
// the largest local signal -- the right choice for a first slice. A
// second-order flux (MUSCL/Roe/HLLC) is a natural next step, not part of
// this class's own scope.
//
// **Time integration: forward Euler**, matching Rusanov's own first order
// in space -- the classic Godunov pairing. stableTimeStep() carries its own
// CFL safety factor internally, the same convention ExplicitTimeStep.hpp
// uses for the incompressible solvers, but this class does not extend that
// shared file: only one solver needs this formula so far, and that file's
// own history (DIVIDA_TECNICA.md 4.1) is what a *duplicated* formula costs
// -- extracting before a second caller exists is the mirror-image mistake.
class CompressibleEulerSolver1D {
public:
    // Chosen once for both ends -- this class has no case yet that needs
    // the two ends to differ.
    enum class BoundaryCondition {
        // Ghost cell mirrors the adjacent real cell's state unchanged --
        // the open-ended treatment a shock tube needs, valid as long as no
        // wave reaches the boundary before the run ends.
        Transmissive,
        // Ghost cell mirrors rho and p, negates velocity -- forces u=0
        // exactly at the face by symmetry of the flux there. This is what
        // makes mass and energy conservation exact unconditionally (see
        // totalMass()/totalEnergy()), not just an alternative boundary.
        Reflecting,
    };

    CompressibleEulerSolver1D(std::size_t n, double length, double gamma, BoundaryCondition boundary);

    // Evaluated once per cell centroid -- same convention as
    // UnstructuredDiffusionSolver::setSourceTerm/setConductivity.
    // `profile(x, density, velocity, pressure)` writes the primitive state
    // at that point.
    void initialize(
        const std::function<void(double x, double& density, double& velocity, double& pressure)>&
            profile);

    // One forward-Euler step of the Rusanov-flux finite-volume update.
    void step(double dt);

    // CFL-limited step: cfl * dx / max_i(|u_i| + c_i). Measured on the Sod
    // shock tube (n=100): stable through cfl=1.05, NaN at cfl=1.10 --
    // consistent with the classic CFL<=1 bound for an explicit Godunov-type
    // scheme, not merely assumed from a textbook number. The default carries
    // roughly 2.5x margin below that measured break, the same kind of
    // conservative gap ExplicitTimeStep.hpp's own 0.3 keeps for the
    // incompressible solvers.
    double stableTimeStep(double cfl = 0.4) const;

    double density(std::size_t cell) const { return rho_.at(cell); }
    double velocity(std::size_t cell) const { return rhoU_.at(cell) / rho_.at(cell); }
    double pressure(std::size_t cell) const;
    double soundSpeed(std::size_t cell) const;
    double cellCenter(std::size_t cell) const;
    double time() const { return time_; }
    std::size_t cellCount() const { return n_; }

    // Domain integrals (sum of rho*dx, rho*u*dx, E*dx) -- what the
    // conservation test checks. Exact (up to floating point) for the
    // whole run with BoundaryCondition::Reflecting: the mass and energy
    // flux at a reflecting wall is exactly zero because u=0 there by
    // construction, for any initial condition, no symmetry required.
    // Not true for Transmissive (an open end generally carries flux) or,
    // even with Reflecting, for momentum (a wall exerts a real net force
    // whenever the two ends' pressures differ).
    double totalMass() const;
    double totalMomentum() const;
    double totalEnergy() const;

private:
    struct State {
        double rho;
        double rhoU;
        double E;
    };

    // Ghost-aware state lookup: i in [0, n) is a real cell, i == -1 and
    // i == n resolve through the boundary condition. The same
    // mirror-by-signed-index shape MixingLengthChannelFlowSolver1D::uAt
    // already uses for its own ghost cells.
    State stateAt(long long i) const;

    static double pressureOf(const State& s, double gamma);
    static double soundSpeedOf(const State& s, double gamma);
    static State flux(const State& s, double gamma);
    static State rusanovFlux(const State& left, const State& right, double gamma);

    std::size_t n_;
    double length_;
    double dx_;
    double gamma_;
    BoundaryCondition boundary_;
    double time_ = 0.0;
    std::vector<double> rho_;
    std::vector<double> rhoU_;
    std::vector<double> E_;
};

} // namespace aether::solver
