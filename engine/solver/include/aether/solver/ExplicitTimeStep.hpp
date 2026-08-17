#pragma once

#include <initializer_list>

namespace aether::solver {

// The explicit stability limit shared by every forward-Euler solver in this
// project, in one place and **with the safety factor already applied**.
//
// **Why this exists, and why the factor is inside rather than at the call
// site.** Both halves of that sentence were the defect recorded as
// DIVIDA_TECNICA.md 4.1.
//
// The formula -- min(von Neumann diffusive limit, CFL convective limit) --
// was copied into seven solvers. That is the situation item 2.1 dealt with
// on the unstructured side, and it had already bitten here: one of the seven
// copies, TransientDiffusionSolver's, was missing the factor of two in
// 1/(2 nu sum 1/h^2) and so returned **twice** the stability limit. Measured
// rather than argued: a 41-cell sine profile stepped at exactly that value
// reaches NaN in 654 steps.
//
// The factor was inside every *caller* instead: forty-nine call sites across
// the tests, the viewer app and the turbulent solvers' primers each wrote
// `0.3 * solver.stableTimeStep()`. A function named stableTimeStep() that
// returns a step which is not stable is a trap, and it was already sprung --
// engine/analysis and engine/persistence both call
// `solver.step(solver.stableTimeStep())` with no factor, having read the
// name and believed it. Those two suites were running at CFL 1.0000, which
// is precisely the fragility Fase 1 diagnosed and this item recorded.
//
// So: one formula, the factor inside it, and the callers stripped of their
// private copies of the margin. The value 0.3 is what forty-nine independent
// call sites had converged on, so moving it inside leaves those cases
// bit-for-bit unchanged -- which is what makes the change checkable.
inline constexpr double kExplicitStabilitySafety = 0.3;

// `spacings` lists the grid spacing along each direction the solver actually
// resolves; a direction with a single cell is simply left out, since it
// carries no second difference and would otherwise dominate the limit.
//
// `velocityScale` is the largest speed the convective term has to carry.
// Pass 0 for a problem with no convection, which drops the CFL term rather
// than dividing by zero.
//
// `effectiveViscosity` is the molecular value plus any eddy viscosity a
// turbulence closure has produced: the diffusive bound is set by what the
// scheme actually diffuses with, not by the fluid's own property.
double explicitStableTimeStep(double effectiveViscosity, double velocityScale,
                               std::initializer_list<double> spacings);

} // namespace aether::solver
