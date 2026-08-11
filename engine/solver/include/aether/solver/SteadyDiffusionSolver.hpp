#pragma once

#include "aether/solver/DiffusionProblem.hpp"

#include <cstddef>
#include <vector>

namespace aether::solver {

// Solves the steady-state Poisson equation nabla^2(phi) = -source (see
// DiffusionProblem for the shared boundary-condition/discretization
// details). source = 0 gives Laplace's equation (e.g. steady heat
// conduction with no heat generation); source != 0 covers anything
// reducible to the same equation form -- including, notably, fully-
// developed steady laminar (Poiseuille) flow between parallel plates,
// where the streamwise momentum equation reduces to mu * d^2u/dy^2 =
// dp/dx, i.e. this same Poisson equation with phi = u and
// source = -(1/mu) * dp/dx. That reuse is why the class and its field are
// named generically (phi/"value"), not "temperature": it is not only a
// heat solver.
//
// Two solve methods are provided: solve() (Gauss-Seidel, simple but slow
// to converge) and solveConjugateGradient() (much faster for this SPD
// system; both solve the identical problem and should agree). This is the
// first working solver Module 5 of the plan calls for; Multigrid/GMRES/
// BiCGSTAB/preconditioners, the full incompressible Navier-Stokes system
// (pressure-velocity coupling, convection), and other physics are later
// work built on top of this same mesh/field layer. See
// TransientDiffusionSolver for the unsteady counterpart of this same
// equation.
class SteadyDiffusionSolver : public DiffusionProblem {
public:
    using DiffusionProblem::DiffusionProblem;

    // Runs up to maxIterations Gauss-Seidel sweeps, stopping early once the
    // largest per-cell change drops below tolerance. Returns the number of
    // iterations actually run.
    std::size_t solve(std::size_t maxIterations = 10000, double tolerance = 1e-9);

    // Conjugate Gradient, matrix-free (each iteration applies the same
    // finite-volume Poisson operator solve() uses, via applyOperator(),
    // rather than building an explicit sparse matrix). Valid because the
    // discretized system is symmetric positive-definite. Converges in far
    // fewer iterations than Gauss-Seidel for the same problem. Stops early
    // once the residual norm drops below tolerance; returns the number of
    // iterations actually run.
    std::size_t solveConjugateGradient(std::size_t maxIterations = 10000, double tolerance = 1e-9);

    // Jacobi-preconditioned Conjugate Gradient: same matrix-free operator,
    // but each iteration also applies M^-1 = 1/diag(A) (the operator's own
    // diagonal -- weightTotal at each free cell, 1 at fixed cells) to the
    // residual before using it as the search direction. The cheapest
    // possible preconditioner (no extra matrix, just one extra vector),
    // and correspondingly modest: it helps most when the diagonal varies a
    // lot between cells (e.g. genuinely anisotropic dx/dy/dz spacing),
    // less on the near-uniform grids this project's tests use so far --
    // see solver_tests.cpp for the measured iteration-count comparison.
    std::size_t solvePreconditionedConjugateGradient(std::size_t maxIterations = 10000,
                                                       double tolerance = 1e-9);

private:
    // Matrix-free application of the discrete Poisson operator: for free
    // cells, weightTotal * x[i] - sum(weight * x[neighbor]); for fixed
    // (Dirichlet) cells, the identity x[i] -- so the linear system's fixed
    // rows simply assert "x equals its prescribed value", which
    // solveConjugateGradient() relies on to keep those entries untouched.
    std::vector<double> applyOperator(const std::vector<double>& x) const;

    // The operator's own diagonal (weightTotal per free cell, 1 for fixed
    // cells) -- the Jacobi preconditioner is just the reciprocal of this.
    std::vector<double> computeDiagonal() const;

    static double dot(const std::vector<double>& a, const std::vector<double>& b);
};

} // namespace aether::solver
