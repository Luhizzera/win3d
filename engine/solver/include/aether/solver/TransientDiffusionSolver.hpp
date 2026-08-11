#pragma once

#include "aether/solver/DiffusionProblem.hpp"

#include <cstddef>

namespace aether::solver {

// Solves the transient (unsteady) diffusion equation
// d(phi)/dt = nabla^2(phi) + source explicitly (forward Euler) on a
// StructuredGrid3D. Boundary conditions and the source-term convention
// are shared with SteadyDiffusionSolver (see DiffusionProblem) - this
// equation's steady state (d(phi)/dt = 0) is exactly
// SteadyDiffusionSolver's nabla^2(phi) = -source, so the two agree at
// t -> infinity for the same boundary conditions and source.
//
// Forward Euler is conditionally stable: step(dt) is only accurate/stable
// for dt no larger than stableTimeStep() (the classic explicit-diffusion
// von Neumann limit, 1/(1/hx^2 + 1/hy^2 + 1/hz^2) over whichever axes
// actually have more than one cell); callers should use a safety margin
// below that (e.g. half), not the exact marginal value. Implicit
// (unconditionally stable) time-stepping is future work.
class TransientDiffusionSolver : public DiffusionProblem {
public:
    using DiffusionProblem::DiffusionProblem;

    // Sets a cell's value directly, regardless of whether it is fixed by
    // setBoundaryValue(). Used to set the initial condition on free cells
    // before stepping; if called on a cell already fixed by
    // setBoundaryValue(), it silently changes that cell's Dirichlet value
    // going forward too, since step() only ever leaves fixed cells alone.
    void setValue(std::size_t i, std::size_t j, std::size_t k, double newValue);

    // The largest forward-Euler step guaranteed von Neumann stable for
    // pure diffusion on this grid.
    double stableTimeStep() const;

    // Advances the field by one explicit forward-Euler step of size dt
    // and advances time() by dt.
    void step(double dt);

    double time() const { return time_; }

private:
    double time_ = 0.0;
};

} // namespace aether::solver
