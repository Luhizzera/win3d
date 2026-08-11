#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 6 (turbulence), third closure: fully-developed turbulent channel
// flow via Menter's k-omega SST (Shear Stress Transport) model -- the
// natural next step after KEpsilonChannelFlowSolver1D, on the same 1D
// channel problem, same exact-force-balance friction velocity, same
// warm-start-from-mixing-length and Patankar-linearization practices.
//
// SST blends two model sets with a wall-distance-dependent function F1:
// Wilcox's k-omega close to the wall (set 1: sigma_k1=0.85, sigma_w1=0.5,
// beta1=0.075) and a k-epsilon-derived k-omega form away from it (set 2:
// sigma_k2=1.0, sigma_w2=0.856, beta2=0.0828), with beta*=0.09 (same value
// as k-epsilon's C_mu) and gamma_i = beta_i/beta* - sigma_wi*kappa^2/sqrt(beta*)
// for each set. A second function F2 gates a Bradshaw shear-stress limiter
// on the eddy viscosity: nu_t = a1*k / max(a1*omega, |dU/dy|*F2), a1=0.31.
// These specific constants are standard, widely-published values (Menter
// 1994/2003) recalled from the literature -- like the log-law's B=5.0
// constant already used by KEpsilonChannelFlowSolver1D, this is *not*
// something derivable from first principles, and is the main place this
// class's correctness cannot be independently cross-checked the way most
// of this project's solvers are. Validation is therefore scoped to what
// *is* self-checkable: the model's own log-law slope and profile
// symmetry (see KEpsilonChannelFlowSolver1D's equivalent caveat).
//
// Governing equations (steady, fully-developed, y-only):
//   0 = P_k - beta* * k * omega   + d/dy[(nu + sigma_k*nu_t)   dk/dy]
//   0 = (gamma/nu_t)*P_k - beta*omega^2 + d/dy[(nu + sigma_w*nu_t) domega/dy]
//       + 2*(1-F1)*sigma_w2/omega * dk/dy * domega/dy   (cross-diffusion)
// with P_k = nu_t*(dU/dy)^2, limited to min(P_k, 10*beta**k*omega) (the
// standard SST production limiter, guards runaway production).
//
// Same wall treatment as KEpsilonChannelFlowSolver1D: equilibrium wall
// functions (mesh does not resolve to the viscous sublayer), k and omega
// fixed at the wall-adjacent cells to k_wall = u_tau^2/sqrt(beta*) and
// omega_wall = u_tau/(kappa*y_wall*sqrt(beta*)) (the latter derived from
// consistency with epsilon = beta**k*omega and KEpsilonChannelFlowSolver1D's
// own eps_wall formula, not a separately memorized value), and u fixed at
// the log-law value (same null-space reasoning as KEpsilonChannelFlowSolver1D:
// a flux condition at both walls leaves momentum with no Dirichlet anchor).
//
// Same Patankar-linearization and warm-start practices as
// KEpsilonChannelFlowSolver1D, applied to two bilinear destruction terms
// instead of one: k's destruction (beta**k*omega, a genuine product of two
// unknowns) is linearized using the previous sweep's omega as the frozen
// factor; omega's destruction (beta*omega^2, the same quadratic-sink shape
// as k-epsilon's epsilon^2/k term) uses the previous sweep's omega as one
// frozen factor. nu_t is under-relaxed for the same feedback-loop-stability
// reason as k-epsilon's nu_t. The velocity field is warm-started from
// MixingLengthChannelFlowSolver1D to avoid the same cold-start production
// deadlock found there (u=0 gives zero production, so k/omega collapse
// before u ever develops).
class KOmegaSSTChannelFlowSolver1D {
public:
    KOmegaSSTChannelFlowSolver1D(std::size_t ny, double height, double kinematicViscosity, double source);

    std::size_t solve(std::size_t maxOuterIterations = 2000, double tolerance = 1e-10);

    double u(std::size_t j) const { return u_[j]; }
    double k(std::size_t j) const { return k_[j]; }
    double omega(std::size_t j) const { return omega_[j]; }
    double eddyViscosity(std::size_t j) const { return nut_[j]; }
    double wallDistance(std::size_t j) const;
    double cellCenterY(std::size_t j) const;

    // Exact friction velocity from the integral momentum balance -- see
    // MixingLengthChannelFlowSolver1D. Independent of k, omega, nu_t.
    double frictionVelocity() const;

private:
    double velocityGradientAt(std::size_t j) const; // interior cells only (1 <= j <= ny_-2)
    void updateBlendingAndCoefficients();

    std::size_t ny_;
    double height_;
    double h_;
    double nu_;
    double source_;
    std::vector<double> u_;
    std::vector<double> k_;
    std::vector<double> omega_;
    std::vector<double> nut_;
    std::vector<double> production_;   // P_k, already limited
    std::vector<double> crossDiffusion_;
    std::vector<double> sigmaK_;
    std::vector<double> sigmaOmega_;
    std::vector<double> beta_;
    std::vector<double> gamma_;
};

} // namespace aether::solver
