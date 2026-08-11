#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 6 (turbulence), second stage: fully-developed turbulent channel
// flow closed with the standard two-equation k-epsilon model, the natural
// next step after MixingLengthChannelFlowSolver1D (same 1D channel
// problem, same exact-force-balance friction velocity, same Picard/Gauss-
// Seidel iterative style) -- but now with transported turbulent kinetic
// energy k and dissipation rate epsilon instead of an algebraic mixing
// length.
//
// Standard model constants: C_mu=0.09, C_eps1=1.44, C_eps2=1.92,
// sigma_k=1.0, sigma_eps=1.3, kappa=0.41 (von Karman).
// nu_t = C_mu * k^2 / epsilon; production P_k = nu_t * (du/dy)^2.
//
// Wall treatment uses equilibrium wall functions, standard for high-
// Reynolds k-epsilon: the mesh does *not* resolve down to the viscous
// sublayer -- the first/last cell centers sit in the log-law region
// (y+ roughly 30-100). k and epsilon there are *fixed* (not solved) to
// the equilibrium values k_wall = u_tau^2/sqrt(C_mu),
// eps_wall = u_tau^3/(kappa*y_wall); u there is *fixed* to the log-law
// value u_tau*((1/kappa)*ln(y+) + B), B=5.0. Fixing u this way (rather
// than a wall-shear flux/Neumann condition) is required, not just
// convenient: a flux condition at *both* walls leaves the momentum
// equation with no Dirichlet anchor anywhere, so u is only determined up
// to an arbitrary additive constant (the same null-space issue the
// pressure solves in TaylorGreenVortexSolver2D/LidDrivenCavitySolver2D
// have, fixed there by pinning one cell) -- caught here because an
// earlier flux-BC version produced an otherwise-sensible but uniformly,
// unphysically negative velocity profile. Using the log law to anchor
// these two cells pulls in the empirical constant B, the one part of this
// project's turbulence validation that *is* a recalled literature value
// (~5.0) rather than derived or measured -- but only to fix the profile's
// absolute level; the log-law *slope* the test checks does not depend on
// having B exactly right.
//
// One simplification specific to validating this class, *not* something a
// general-purpose solver could do: u_tau is obtained from the exact
// integral momentum balance (sqrt(source*height/2), independent of the
// turbulence closure -- see MixingLengthChannelFlowSolver1D), rather than
// iteratively inverting the log law against a resolved near-wall velocity
// (what real wall-function implementations must do for a general
// geometry where the wall shear stress isn't known in advance).
//
// All source/production/destruction terms are refreshed once per outer
// Picard iteration and held frozen through that iteration's inner Gauss-
// Seidel sweeps (same style as the mixing-length solver, just applied to
// three coupled fields -- u, k, epsilon -- instead of one). nu_t itself is
// under-relaxed (factor 0.3): nu_t=Cmu*k^2/epsilon amplifies noise in k
// and epsilon, and without relaxation this feedback loop oscillated
// rather than settling during development. The epsilon equation's
// quadratic destruction term is linearized implicitly (Patankar's
// positivity rule -- sink terms in the denominator, not the numerator):
// treating it fully explicitly caused runaway collapse of k and epsilon
// to zero during development, a documented, common instability for this
// kind of stiff quadratic sink rather than something specific to this
// implementation.
//
// The velocity field is warm-started from MixingLengthChannelFlowSolver1D
// run on the same mesh/parameters, rather than from rest. Starting from
// u=0 gives zero production (P_k = nu_t*(du/dy)^2 = 0 with no velocity
// gradient anywhere), and with nothing to balance the destruction term k
// collapses toward zero within the first few outer iterations, before the
// velocity field ever gets a chance to develop -- a real chicken-and-egg
// problem this class hit during development. The mixing-length model
// already gives a reasonable, cheap, validated-elsewhere velocity profile
// to break that deadlock.
class KEpsilonChannelFlowSolver1D {
public:
    KEpsilonChannelFlowSolver1D(std::size_t ny, double height, double kinematicViscosity, double source);

    std::size_t solve(std::size_t maxOuterIterations = 2000, double tolerance = 1e-10);

    double u(std::size_t j) const { return u_[j]; }
    double k(std::size_t j) const { return k_[j]; }
    double epsilon(std::size_t j) const { return epsilon_[j]; }
    double eddyViscosity(std::size_t j) const { return nut_[j]; }
    double wallDistance(std::size_t j) const;
    double cellCenterY(std::size_t j) const;

    // Exact friction velocity from the integral momentum balance -- see
    // the class comment. Independent of k, epsilon, nu_t.
    double frictionVelocity() const;

private:
    double velocityGradientAt(std::size_t j) const; // interior cells only (1 <= j <= ny_-2)
    void updateEddyViscosityAndProduction();

    std::size_t ny_;
    double height_;
    double h_;
    double nu_;
    double source_;
    std::vector<double> u_;
    std::vector<double> k_;
    std::vector<double> epsilon_;
    std::vector<double> nut_;
    std::vector<double> production_;
};

} // namespace aether::solver
