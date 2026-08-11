#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 6, extending the two-equation k-epsilon closure (previously only
// coupled to the 1D fully-developed channel, KEpsilonChannelFlowSolver1D)
// to a genuinely 2D convecting/recirculating flow -- the lid-driven cavity,
// same domain MixingLengthLidDrivenCavitySolver2D already validated a
// simpler (algebraic) closure against. This is a materially harder step
// than that one: k and epsilon are now transported fields (their own
// advection-diffusion-production-destruction equations, explicitly time-
// marched alongside momentum every step, not solved to convergence by an
// outer Picard loop the way the 1D channel version was), and a cavity has
// four walls with no single "streamwise" direction, unlike a channel's two
// parallel walls.
//
// **Wall treatment is the central design choice, and deliberately not
// equilibrium wall functions** (what KEpsilonChannelFlowSolver1D uses):
// equilibrium wall functions need the local friction velocity u_tau, which
// the 1D channel got from an *exact* integral force balance specific to
// fully-developed flow (sqrt(source*height/2)) -- no such exact relation
// exists for a general 2D recirculating flow, and estimating u_tau locally
// from a resolved near-wall gradient on a mesh not specifically refined for
// it would be exactly the kind of unvalidated approximation this project
// avoids. Instead this class uses the low-Reynolds-number *asymptotic*
// near-wall relation for epsilon, derived from the Taylor-series fact that
// k ~ y^2 approaching any solid wall (true generally, not just for
// log-law/equilibrium flows): epsilon_wall = 2*nu*k_nearWallCell/y_half^2,
// where y_half = h/2 is the exact, known distance from a cell center to
// the wall under this class's ghost-mirror convention. k itself is fixed
// to exactly 0 at every wall (turbulent kinetic energy vanishes at a
// no-slip wall by definition). Cmu is *not* damped near the wall (no f_mu
// low-Re correction, e.g. Launder-Sharma) -- a further simplification,
// since that would pull in more recalled empirical functions; accuracy
// very close to the wall on this class's modest mesh resolution is not
// claimed.
//
// Same ghost-mirror Dirichlet treatment, tapered lid, Chorin projection,
// and warm-start practice as MixingLengthLidDrivenCavitySolver2D (velocity
// is warm-started from that simpler closure, for the same cold-start-
// deadlock reason KEpsilonChannelFlowSolver1D warm-starts from
// MixingLengthChannelFlowSolver1D: u=0 gives zero production, so k and
// epsilon would collapse before the velocity field ever develops).
class KEpsilonLidDrivenCavitySolver2D {
public:
    KEpsilonLidDrivenCavitySolver2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                                     double viscosity, double lidVelocity);

    double stableTimeStep() const;
    void step(double dt);

    double u(std::size_t i, std::size_t j) const;
    double v(std::size_t i, std::size_t j) const;
    double pressure(std::size_t i, std::size_t j) const;
    double k(std::size_t i, std::size_t j) const;
    double epsilon(std::size_t i, std::size_t j) const;
    double eddyViscosity(std::size_t i, std::size_t j) const;
    double time() const { return time_; }

    double maxDivergence() const;

private:
    std::size_t index(std::size_t i, std::size_t j) const { return i + j * nx_; }

    double lidVelocityAt(std::size_t i) const;

    double dirichletAt(const std::vector<double>& field, std::size_t i, std::size_t j, int di, int dj,
                        double topWallValue) const;
    double neumannAt(const std::vector<double>& field, std::size_t i, std::size_t j, int di, int dj) const;
    // epsilon's ghost mirrors a *computed* wall value (2*nu*k/y_half^2,
    // using the adjacent interior cell's own k and the known half-cell
    // wall distance) -- see the class comment.
    double epsilonGhostAt(std::size_t i, std::size_t j, int di, int dj) const;
    // Eddy viscosity at a neighbor; exactly 0.0 across any wall crossing
    // (nu_t vanishes at a solid wall).
    double nutAt(std::size_t i, std::size_t j, int di, int dj) const;

    void updateEddyViscosityAndProduction();

    std::vector<double> applyLaplacian(const std::vector<double>& x) const;
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    void projectToDivergenceFree(std::vector<double>& uStar, std::vector<double>& vStar, double dt);

    std::size_t nx_;
    std::size_t ny_;
    double lengthX_;
    double lengthY_;
    double dx_;
    double dy_;
    double viscosity_;
    double lidVelocity_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> p_;
    std::vector<double> k_;
    std::vector<double> epsilon_;
    std::vector<double> nut_;
    std::vector<double> production_;
    double time_ = 0.0;
};

} // namespace aether::solver
