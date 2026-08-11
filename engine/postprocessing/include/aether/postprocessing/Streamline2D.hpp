#pragma once

#include "aether/core/Vector3.hpp"

#include <cstddef>
#include <vector>

namespace aether::postprocessing {

// Module 7: real streamline integration through a 2D velocity field --
// distinct from what apps/cavity_viewer draws today (arrows at cell
// centers, not integrated paths). Traces a curve tangent everywhere to the
// (interpolated) velocity field via 4th-order Runge-Kutta, the standard
// choice for streamline integration (much lower numerical drift per step
// than Euler, without needing an implicit solve).
//
// The velocity field is given as two nx*ny arrays, cell-centered (matching
// the convention every collocated 2D Navier-Stokes solver in this project
// already uses: u(i,j) at ((i+0.5)*dx, (j+0.5)*dy)), sampled at arbitrary
// (non-grid-aligned) points via bilinear interpolation between the 4
// nearest cell centers.
class Streamline2D {
public:
    // periodic=true wraps sample points modulo the domain (matching
    // TaylorGreenVortexSolver2D's convention); periodic=false clamps to
    // the domain edge (matching LidDrivenCavitySolver2D's solid walls) --
    // trace() also stops early in that case if the path would leave the
    // domain, rather than clamping the streamline itself to the boundary.
    Streamline2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                 const std::vector<double>& u, const std::vector<double>& v, bool periodic);

    // Velocity at an arbitrary point via bilinear interpolation of the 4
    // nearest cell centers.
    core::Vector3 velocityAt(double x, double y) const;

    // Integrates forward from (x0, y0) via RK4 with fixed step size
    // stepSize (physical units, not cell count), for up to maxSteps
    // points. Stops early (returning fewer points) if velocity magnitude
    // drops below 1e-12 (a stagnation point) or, for a non-periodic
    // domain, if the path would leave [0,lengthX]x[0,lengthY].
    std::vector<core::Vector3> trace(double x0, double y0, double stepSize, std::size_t maxSteps) const;

private:
    std::size_t nx_;
    std::size_t ny_;
    double lengthX_;
    double lengthY_;
    double dx_;
    double dy_;
    bool periodic_;
    std::vector<double> u_;
    std::vector<double> v_;
};

} // namespace aether::postprocessing
