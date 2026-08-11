#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Solves the 2D incompressible Navier-Stokes equations (constant density,
// kinematic viscosity nu) on a doubly-periodic collocated grid via Chorin's
// projection method:
//   1. Explicit predictor: advance u, v ignoring pressure, using central
//      differences for both convection and diffusion.
//   2. Pressure projection: solve nabla^2(p) = div(u*)/dt (Conjugate
//      Gradient, periodic 5-point Laplacian, one cell pinned to remove the
//      null space every periodic/pure-Neumann Poisson problem has), then
//      subtract dt*grad(p) from the predicted velocity to make it
//      divergence-free again.
//
// Deliberately scoped to a doubly-periodic domain (no solid walls) so it
// can be validated against the Taylor-Green vortex, an EXACT closed-form
// solution of the full nonlinear incompressible Navier-Stokes equations
// (unlike the Poiseuille flow SteadyDiffusionSolver validates against,
// which is a fully-developed special case where convection drops out).
// This is the project's first real Navier-Stokes solver with convection
// and pressure-velocity coupling. Solid-wall boundary conditions (e.g. for
// the classic lid-driven cavity benchmark) need a staggered grid and
// ghost-cell wall treatment not implemented here -- future work.
//
// Uses a collocated (not staggered) arrangement -- u, v, p all cell-
// centered -- which is simpler to implement than the classic staggered
// SIMPLE layout but can checkerboard on coarser/rougher problems without
// Rhie-Chow interpolation (not implemented). Acceptable for this smooth,
// well-resolved periodic test case; a known limitation for general use.
class TaylorGreenVortexSolver2D {
public:
    TaylorGreenVortexSolver2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                               double viscosity);

    // Sets the velocity at a cell (used to set the initial condition).
    void setVelocity(std::size_t i, std::size_t j, double u, double v);

    // A conservative explicit-stability estimate combining the diffusive
    // von Neumann limit with a convective CFL limit based on the given
    // characteristic velocity scale. Callers should use a margin below
    // this (e.g. half), not the exact marginal value.
    double stableTimeStep(double velocityScale) const;

    // Advances one step (predictor + pressure projection) and advances
    // time() by dt.
    void step(double dt);

    double u(std::size_t i, std::size_t j) const;
    double v(std::size_t i, std::size_t j) const;
    double pressure(std::size_t i, std::size_t j) const;
    double time() const { return time_; }

    // Maximum absolute discrete divergence over all cells for the
    // *current* velocity field. Should be near zero after step(), since
    // removing divergence is exactly what the pressure projection does --
    // a direct, self-contained correctness check that needs no external
    // reference data.
    double maxDivergence() const;

private:
    std::size_t index(std::size_t i, std::size_t j) const { return i + j * nx_; }
    std::size_t wrap(long long i, std::size_t n) const;

    std::vector<double> applyPeriodicLaplacian(const std::vector<double>& x) const;
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    void projectToDivergenceFree(std::vector<double>& uStar, std::vector<double>& vStar, double dt);

    std::size_t nx_;
    std::size_t ny_;
    double dx_;
    double dy_;
    double viscosity_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> p_;
    double time_ = 0.0;
};

} // namespace aether::solver
