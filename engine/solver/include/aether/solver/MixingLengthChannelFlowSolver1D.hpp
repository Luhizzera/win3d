#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 6 (turbulence), first pass: fully-developed turbulent channel
// flow closed with Prandtl's mixing-length model -- the historically first
// and simplest RANS turbulence closure (an algebraic/"zero-equation"
// model, no extra transport PDEs), used here deliberately before
// attempting full two-equation k-epsilon. Getting the wall treatment and
// log-law validation right with this simpler model first, on the exact
// same 1D fully-developed channel problem SteadyDiffusionSolver already
// validates laminar (Poiseuille) flow against, is lower-risk than adding
// k-epsilon's two coupled nonlinear transport equations and near-wall
// damping/wall-function machinery in one step. k-epsilon is the natural
// next stage once this is proven.
//
// Governing equation (steady, fully-developed, so convection vanishes and
// only y varies): 0 = source + d/dy[(nu + nu_t) du/dy], with no-slip
// (u=0) at both walls y=0 and y=height. source stands for -(1/rho)(dp/dx),
// the same convention as SteadyDiffusionSolver's Poiseuille validation.
//
// Turbulence closure: nu_t(y) = l_m(y)^2 * |du/dy|, with mixing length
// l_m(y) = min(kappa * d(y), 0.09 * height/2) -- kappa=0.41 (von Karman
// constant) times distance to the nearest wall, capped at the classic
// Escudier asymptotic value far from the wall. Solved by Picard (fixed-
// point) iteration: freeze nu_t, solve the resulting linear variable-
// coefficient diffusion equation by Gauss-Seidel, recompute nu_t from the
// new velocity gradient, repeat.
//
// No van Driest near-wall damping and no wall functions are implemented,
// so the viscous sublayer (y+ < ~5) is not expected to be accurate; only
// the log-law region (y+ roughly 30-100) is validated against, where the
// mixing-length ansatz's log-law slope 1/kappa is close to guaranteed by
// construction rather than an independent empirical fact being recalled
// from memory.
class MixingLengthChannelFlowSolver1D {
public:
    MixingLengthChannelFlowSolver1D(std::size_t ny, double height, double kinematicViscosity,
                                     double source);

    // Picard iteration outer loop (each iteration: update nu_t, then a
    // fixed number of Gauss-Seidel sweeps on u with nu_t frozen). Stops
    // early once the largest per-cell change in u between outer
    // iterations drops below tolerance. Returns the number of outer
    // iterations actually run.
    std::size_t solve(std::size_t maxOuterIterations = 500, double tolerance = 1e-10);

    double u(std::size_t j) const { return u_[j]; }
    double eddyViscosity(std::size_t j) const { return nut_[j]; }
    double wallDistance(std::size_t j) const;
    double cellCenterY(std::size_t j) const;

    // Friction velocity u_tau = sqrt(tau_wall / rho), computed from the
    // *exact* integral momentum balance for fully-developed channel flow
    // (integrating 0 = source + d/dy[(nu+nu_t) du/dy] from the wall to the
    // centerline, using du/dy=0 at the centerline by symmetry and nu_t=0
    // at the wall by definition of the mixing-length model): this gives
    // tau_wall/rho = source * (height/2) exactly, regardless of the
    // turbulence closure's details. Deliberately not estimated from a
    // finite-difference wall-gradient of the discrete solution: an earlier
    // version did that using pure molecular viscosity, inconsistent with
    // solve()'s actual wall-face diffusivity (which blends in the
    // adjacent cell's nu_t) -- a ~38% error caught by comparing against
    // this exact formula, which has no such discretization ambiguity.
    double frictionVelocity() const;

private:
    double uAt(long long j) const; // Dirichlet-mirrored across either wall
    double gammaAt(long long j) const; // nu + nu_t, mirrored across either wall
    double velocityGradientAt(std::size_t j) const;
    void updateEddyViscosity();

    std::size_t ny_;
    double height_;
    double h_;
    double nu_;
    double source_;
    std::vector<double> u_;
    std::vector<double> nut_;
};

} // namespace aether::solver
