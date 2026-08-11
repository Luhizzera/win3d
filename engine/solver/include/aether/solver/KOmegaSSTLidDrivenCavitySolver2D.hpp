#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 6: k-omega SST (previously only coupled to the 1D fully-developed
// channel, KOmegaSSTChannelFlowSolver1D) extended to the same real 2D
// convecting lid-driven cavity KEpsilonLidDrivenCavitySolver2D already
// validated the simpler two-equation closure against. Mirrors that class's
// architecture closely (ghost-mirror walls, explicit transport stepping,
// warm-start from mixing length, Patankar-linearized destruction terms)
// while swapping in SST's blended F1/F2 formulation.
//
// **Wall treatment for omega is self-derived, not the standard Wilcox
// low-Re formula recalled from memory.** KEpsilonLidDrivenCavitySolver2D
// derives its epsilon wall condition from the general near-wall Taylor-
// series fact k ~ y^2 (giving epsilon_wall = 2*nu*k/y_half^2). Substituting
// that into the standard, exact k-omega/k-epsilon conversion identity
// epsilon = beta* * k * omega (valid throughout the flow) and taking the
// y -> 0 limit (where k ~ a*y^2 makes the k-dependence cancel) gives
// omega_wall = 2*nu / (beta* * y_half^2) -- finite in form, independent of
// the local k value, and consistent with the already-validated epsilon
// relation. This is *not* claimed to exactly match Wilcox's more careful
// matched-asymptotic near-wall formula (commonly cited with a factor of 6
// rather than 2, from directly analyzing the omega equation's own near-
// wall balance rather than converting from epsilon) -- that distinction is
// a genuine simplification of this class, documented rather than hidden.
// Every other SST model constant here (sigma_k1/2, sigma_w1/2, beta1/2,
// a1, kappa, beta*) is the same recalled-literature value already flagged
// in KOmegaSSTChannelFlowSolver1D's own comment.
class KOmegaSSTLidDrivenCavitySolver2D {
public:
    KOmegaSSTLidDrivenCavitySolver2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                                      double viscosity, double lidVelocity);

    double stableTimeStep() const;
    void step(double dt);

    double u(std::size_t i, std::size_t j) const;
    double v(std::size_t i, std::size_t j) const;
    double pressure(std::size_t i, std::size_t j) const;
    double k(std::size_t i, std::size_t j) const;
    double omega(std::size_t i, std::size_t j) const;
    double eddyViscosity(std::size_t i, std::size_t j) const;
    double time() const { return time_; }

    double maxDivergence() const;

private:
    std::size_t index(std::size_t i, std::size_t j) const { return i + j * nx_; }

    double lidVelocityAt(std::size_t i) const;
    double wallDistanceAt(std::size_t i, std::size_t j) const;

    double dirichletAt(const std::vector<double>& field, std::size_t i, std::size_t j, int di, int dj,
                        double topWallValue) const;
    double neumannAt(const std::vector<double>& field, std::size_t i, std::size_t j, int di, int dj) const;
    // omega's ghost mirrors a *computed* wall value (2*nu/(beta**y_half^2),
    // independent of local k -- see the class comment).
    double omegaGhostAt(std::size_t i, std::size_t j, int di, int dj) const;
    // Eddy viscosity at a neighbor; exactly 0.0 across any wall crossing.
    double nutAt(std::size_t i, std::size_t j, int di, int dj) const;
    // sigma_k(or omega) * nu_t at a neighbor, for the transport equations'
    // own diffusion coefficient; exactly 0.0 across any wall crossing
    // (nu_t vanishes there, whatever the locally-blended sigma would be).
    double sigmaKNutAt(std::size_t i, std::size_t j, int di, int dj) const;
    double sigmaOmegaNutAt(std::size_t i, std::size_t j, int di, int dj) const;

    void updateBlendingAndCoefficients();

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
    std::vector<double> omega_;
    std::vector<double> nut_;
    std::vector<double> production_;
    std::vector<double> crossDiffusion_;
    std::vector<double> sigmaK_;
    std::vector<double> sigmaOmega_;
    std::vector<double> beta_;
    std::vector<double> gamma_;
    double time_ = 0.0;
};

} // namespace aether::solver
