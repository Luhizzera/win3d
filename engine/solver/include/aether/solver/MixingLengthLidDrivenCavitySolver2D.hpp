#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 6, extended to a genuinely 2D flow: the lid-driven cavity
// (LidDrivenCavitySolver2D) with a Prandtl mixing-length eddy-viscosity
// closure, evaluated pointwise from the *local* 2D flow field instead of
// the 1D fully-developed channel's single wall-normal profile. Every
// turbulence closure so far (MixingLengthChannelFlowSolver1D,
// KEpsilonChannelFlowSolver1D, KOmegaSSTChannelFlowSolver1D) only ever saw
// a 1D problem where convection vanishes by construction (fully developed
// flow); this is the first closure coupled to a real convecting,
// recirculating velocity field.
//
// nu_t(x,y) = l_m(x,y)^2 * S(x,y), with:
//   l_m = min(kappa * d_wall, 0.09 * min(lengthX, lengthY)/2) -- same
//     von-Karman-near-wall / Escudier-cap formula as the 1D model, now
//     using distance to the *nearest* of the cavity's four solid walls
//     (d_wall = min(x, lengthX-x, y, lengthY-y)) instead of a single
//     wall-normal coordinate.
//   S = sqrt(2*(du/dx)^2 + 2*(dv/dy)^2 + (du/dy+dv/dx)^2) -- the full 2D
//     strain-rate magnitude (the 1D channel's |du/dy| is the special case
//     of this where du/dx=dv/dy=dv/dx=0), needed because a cavity has
//     real normal strain (du/dx, dv/dy) in addition to shear.
//
// nu_t is folded into the momentum predictor's diffusion term as
// (nu + nu_t) at each face, face-averaged and applied independently per
// velocity component -- i.e. d/dx[(nu+nu_t) du/dx] + d/dy[(nu+nu_t) du/dy]
// for u, and the same pattern for v. This is the standard simplified
// ("quasi-laminar-form") treatment of a spatially varying eddy viscosity:
// it omits the extra cross-derivative terms a fully rigorous variable-
// viscosity stress-tensor divergence would include (e.g. d/dx[(nu+nu_t)
// dv/dx] contributing to the u-equation). Common in introductory/engineering
// RANS codes and adequate for a first 2D turbulence-coupled solver; a full
// tensor treatment is future work if needed. nu_t is exactly zero at all
// four walls by construction (mixing length vanishes there), matching the
// model's physical definition.
//
// Same solid-wall ghost-cell treatment, tapered lid, Chorin projection
// (unaffected by nu_t -- only the momentum diffusion term uses it, not the
// pressure Poisson operator), and same "no literature benchmark" validation
// philosophy as LidDrivenCavitySolver2D: divergence stays bounded, the
// primary-vortex topology still holds, and nu_t is checked to vanish at
// the walls and be non-negative -- not a claim of resolving fully
// turbulent Re (that would need a much finer mesh and longer transient
// than this class is validated at).
class MixingLengthLidDrivenCavitySolver2D {
public:
    MixingLengthLidDrivenCavitySolver2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                                         double viscosity, double lidVelocity);

    // Conservative estimate based on molecular viscosity and the current
    // (possibly still-zero, if called before any step()) eddy viscosity.
    // Turbulent diffusion grows as the flow develops, so callers should
    // keep a healthy safety margin (e.g. 0.3x, as used throughout this
    // project) rather than relying on this being exact for all time.
    double stableTimeStep() const;

    void step(double dt);

    double u(std::size_t i, std::size_t j) const;
    double v(std::size_t i, std::size_t j) const;
    double pressure(std::size_t i, std::size_t j) const;
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
    // Eddy viscosity at (i+di, j+dj); exactly 0.0 if that neighbor crosses
    // a wall (nu_t vanishes at a solid wall by definition), not mirrored.
    double nutAt(std::size_t i, std::size_t j, int di, int dj) const;

    void updateEddyViscosity();

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
    std::vector<double> nut_;
    double time_ = 0.0;
};

} // namespace aether::solver
