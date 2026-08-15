#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::optimization {

// A black-box scalar cost to minimize over a real parameter vector -- no
// gradient required, so the caller is free to plug in something as
// non-differentiable-in-practice as "build a solver with these parameters,
// run it, read off a derived quantity".
using ObjectiveFunction = std::function<double(const std::vector<double>&)>;

struct OptimizationResult {
    std::vector<double> parameters;
    double value;
    std::size_t iterations;
    bool converged; // false means maxIterations was reached first
};

// Nelder-Mead downhill simplex (Nelder & Mead, 1965): derivative-free
// minimization of an n-dimensional scalar objective. Chosen as Module 12's
// first optimizer because it treats the objective as a pure black box --
// exactly what's needed to optimize *over a CFD solver's own output* (run
// the solver, read off some derived quantity, that is the cost) without
// ever needing the solver's derivative with respect to its parameters.
//
// Standard algorithm: maintain n+1 vertices forming a simplex in parameter
// space; each iteration, replace the worst vertex by reflecting it through
// the centroid of the rest, expanding further if that reflection was a
// large improvement, or contracting toward the centroid if it wasn't an
// improvement at all -- shrinking the whole simplex toward its best vertex
// only when neither reflection nor contraction helps.
class NelderMead {
public:
    explicit NelderMead(std::size_t maxIterations = 1000, double tolerance = 1e-10);

    // `initialStepSize` sizes the initial simplex along each axis from
    // `initialGuess`: too small and it can stall in a shallow local dip
    // near the start, too large wastes early iterations exploring -- the
    // same problem-dependent tradeoff any other iterative method's
    // starting parameters have in this project.
    OptimizationResult minimize(const ObjectiveFunction& objective, const std::vector<double>& initialGuess,
                                 double initialStepSize = 0.1) const;

private:
    std::size_t maxIterations_;
    double tolerance_;
};

} // namespace aether::optimization
