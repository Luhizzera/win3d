#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Solves the 2D incompressible Navier-Stokes equations on a square cavity
// with solid walls: no-slip (u=v=0) on the left, right and bottom walls,
// and a moving lid on top -- the classic lid-driven cavity benchmark. Same
// numerical method as TaylorGreenVortexSolver2D (collocated grid, Chorin
// projection: explicit convection+diffusion predictor, then a Conjugate-
// Gradient pressure projection), but with solid-wall ghost-cell boundary
// treatment instead of periodic wraparound: velocity ghost values are
// mirrored to enforce the exact wall value (Dirichlet, e.g.
// ghost = 2*wallValue - interior), pressure ghost values simply copy the
// interior value (Neumann/zero-gradient, the standard assumption at an
// impermeable wall with no-slip velocity).
//
// The lid's tangential velocity is *regularized* to taper smoothly to zero
// at its two top corners (lidVelocity * sin^2(pi*x/lengthX)) rather than
// holding a constant value with a step discontinuity where it meets the
// stationary side walls. A discontinuous lid is a genuinely singular
// boundary condition for the continuous Navier-Stokes equations (pressure
// blows up at those two corners); measured directly while validating this
// class, that singularity dominated the discrete divergence (an interior
// cell one step away from a top corner had ~20x the divergence of any
// other interior cell) rather than it being spread out as ordinary
// discretization error. Tapering the lid is standard, documented practice
// for exactly this reason, not an ad hoc fix.
//
// Deliberately validated *without* comparing to literature benchmark
// tables (e.g. Ghia et al. 1982): those would have to be recalled from
// memory here, a real accuracy risk. Instead validated against properties
// derivable from first principles: the exact rest state when the lid
// isn't moving, mass conservation (bounded discrete divergence, as with
// TaylorGreenVortexSolver2D), and the flow-reversal that mass conservation
// in a closed box necessarily forces (see solver_tests.cpp for details).
//
// Same collocated-grid limitation as TaylorGreenVortexSolver2D applies:
// checkerboard-prone without Rhie-Chow interpolation (not implemented).
class LidDrivenCavitySolver2D {
public:
    LidDrivenCavitySolver2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                             double viscosity, double lidVelocity);

    double stableTimeStep() const;

    void step(double dt);

    double u(std::size_t i, std::size_t j) const;
    double v(std::size_t i, std::size_t j) const;
    double pressure(std::size_t i, std::size_t j) const;
    double time() const { return time_; }

    double maxDivergence() const;

private:
    std::size_t index(std::size_t i, std::size_t j) const { return i + j * nx_; }

    // The regularized lid speed at cell-column i (see the class comment).
    double lidVelocityAt(std::size_t i) const;

    // Value of `field` at (i+di, j+dj), Dirichlet-mirrored across whichever
    // wall that neighbor would cross (topWallValue is the prescribed value
    // if that wall is the top/lid; every other wall is 0 for both velocity
    // components).
    double dirichletAt(const std::vector<double>& field, std::size_t i, std::size_t j, int di, int dj,
                        double topWallValue) const;

    // Value of `field` at (i+di, j+dj), Neumann-mirrored (ghost equals the
    // boundary cell itself, i.e. zero gradient) across any wall crossed.
    double neumannAt(const std::vector<double>& field, std::size_t i, std::size_t j, int di, int dj) const;

    std::vector<double> applyLaplacian(const std::vector<double>& x) const;
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    void projectToDivergenceFree(std::vector<double>& uStar, std::vector<double>& vStar, double dt);

    std::size_t nx_;
    std::size_t ny_;
    double dx_;
    double dy_;
    double viscosity_;
    double lidVelocity_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> p_;
    double time_ = 0.0;
};

} // namespace aether::solver
