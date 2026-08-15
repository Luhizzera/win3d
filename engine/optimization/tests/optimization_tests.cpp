#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/optimization/NelderMead.hpp"
#include "aether/solver/SteadyDiffusionSolver.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>

using aether::core::Vector3;
using aether::mesh::StructuredGrid3D;
using aether::optimization::NelderMead;
using aether::optimization::OptimizationResult;
using aether::solver::SteadyDiffusionSolver;

namespace {

// Simplest possible exact case: a paraboloid with a known, hand-computable
// global minimum at (3, -2), value 0 -- no literature recall involved, just
// arithmetic anyone can verify (d/dx[(x-3)^2+(y+2)^2] = 2(x-3), zero at
// x=3, symmetric argument for y).
void testNelderMeadMatchesKnownQuadraticMinimum() {
    NelderMead optimizer(500, 1e-14);
    const auto objective = [](const std::vector<double>& p) {
        const double dx = p[0] - 3.0;
        const double dy = p[1] + 2.0;
        return dx * dx + dy * dy;
    };
    const OptimizationResult result = optimizer.minimize(objective, {0.0, 0.0}, 1.0);

    std::printf("  [aether_optimization_tests] quadratic: x=%.10f y=%.10f value=%.3e iters=%zu converged=%d\n",
                result.parameters[0], result.parameters[1], result.value, result.iterations, result.converged);

    AETHER_CHECK(result.converged);
    AETHER_CHECK(std::fabs(result.parameters[0] - 3.0) < 1e-4);
    AETHER_CHECK(std::fabs(result.parameters[1] - (-2.0)) < 1e-4);
    AETHER_CHECK(result.value < 1e-8);
}

// The standard derivative-free-optimizer benchmark: Rosenbrock's function
// f(x,y) = 100*(y-x^2)^2 + (1-x)^2. Its minimum is derivable directly, not
// recalled from a table: both partial derivatives vanish only where
// y = x^2 (from df/dy) and, substituting, (1-x) = 0 (from df/dx), i.e.
// exactly at (1, 1), where f = 0. Famous for a long curved, nearly-flat
// valley that makes derivative-free methods converge slowly along it --
// so this is a genuinely harder case than the quadratic above, not a
// repeat of it.
void testNelderMeadMatchesRosenbrockMinimum() {
    NelderMead optimizer(5000, 1e-15);
    const auto objective = [](const std::vector<double>& p) {
        const double x = p[0];
        const double y = p[1];
        const double a = y - x * x;
        const double b = 1.0 - x;
        return 100.0 * a * a + b * b;
    };
    const OptimizationResult result = optimizer.minimize(objective, {-1.2, 1.0}, 0.5);

    std::printf("  [aether_optimization_tests] rosenbrock: x=%.10f y=%.10f value=%.3e iters=%zu converged=%d\n",
                result.parameters[0], result.parameters[1], result.value, result.iterations, result.converged);

    AETHER_CHECK(std::fabs(result.parameters[0] - 1.0) < 1e-3);
    AETHER_CHECK(std::fabs(result.parameters[1] - 1.0) < 1e-3);
}

// The point of Module 12's first pass: use the optimizer to solve a real
// inverse problem over an already-validated solver, not just a synthetic
// test function. Same 1D "Poiseuille via SteadyDiffusionSolver" setup as
// solver_tests.cpp's testPlanePoiseuilleProfile() (both x-boundaries fixed
// to 0, source term standing in for -(1/mu)*dp/dx): that test already
// proved the exact discrete solution is
// value(i) = (source * h^2 / 2) * i * (nx-1-i). Here the question is run
// backwards -- given a *target* value at one cell, find the source term
// that produces it -- and the optimizer's answer is checked against the
// exact closed-form inverse of that same formula, not just against "the
// residual got small".
void testNelderMeadFindsSourceTermForTargetPoiseuilleValue() {
    const std::size_t nx = 30;
    const double h = 0.1;
    const std::size_t targetCell = 15;

    // value(targetCell) = source * (h*h/2) * targetCell * (nx-1-targetCell)
    //                    = source * 0.005 * 15 * 14 = source * 1.05
    const double coefficient =
        (h * h / 2.0) * static_cast<double>(targetCell) * static_cast<double>(nx - 1 - targetCell);
    const double trueSource = 2.0;
    const double targetValue = trueSource * coefficient;

    const auto objective = [&](const std::vector<double>& p) {
        StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(h * static_cast<double>(nx), 0.05, 0.05), nx, 1, 1);
        SteadyDiffusionSolver solver(grid);
        solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
        solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 0.0);
        solver.setSourceTerm(p[0]);
        solver.solveConjugateGradient(2000, 1e-12);

        const double diff = solver.value(targetCell, 0, 0) - targetValue;
        return diff * diff;
    };

    NelderMead optimizer(200, 1e-20);
    const OptimizationResult result = optimizer.minimize(objective, {1.0}, 0.5);

    std::printf("  [aether_optimization_tests] poiseuille inverse: source=%.10f (true=%.10f) value=%.3e "
                "iters=%zu converged=%d\n",
                result.parameters[0], trueSource, result.value, result.iterations, result.converged);

    AETHER_CHECK(std::fabs(result.parameters[0] - trueSource) < 1e-4);
}

} // namespace

int main() {
    testNelderMeadMatchesKnownQuadraticMinimum();
    testNelderMeadMatchesRosenbrockMinimum();
    testNelderMeadFindsSourceTermForTargetPoiseuilleValue();
    std::printf("aether_optimization_tests: OK\n");
    return 0;
}
