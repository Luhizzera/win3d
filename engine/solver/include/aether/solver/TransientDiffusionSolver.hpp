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
// Forward Euler is conditionally stable: step(dt) is only stable for dt no
// larger than stableTimeStep(), which returns the von Neumann limit
// 1/(2 * sum over resolved axes of 1/h^2) **with a safety factor already
// applied** -- so a caller may use the value as returned, which is the whole
// point of the name.
//
// It did not always mean that. This class used to return 1/(sum 1/h^2),
// missing the factor of two, i.e. exactly twice the stability limit, with a
// comment telling callers to halve it themselves. Measured rather than
// argued when that was found (DIVIDA_TECNICA.md 4.1): a 41-cell sine profile
// stepped at the value it returned reached NaN in 654 steps. The limit now
// lives in explicitStableTimeStep(), shared with every other explicit solver
// here, and the margin is measurable -- the same profile survives up to 3.33x
// the returned value and diverges at 3.40x, which is the 0.3 factor exactly.
//
// Implicit (unconditionally stable) time-stepping is future work.
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
