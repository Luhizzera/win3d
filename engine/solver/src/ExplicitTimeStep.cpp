#include "aether/solver/ExplicitTimeStep.hpp"

#include <algorithm>
#include <limits>

namespace aether::solver {

double explicitStableTimeStep(double effectiveViscosity, double velocityScale,
                               std::initializer_list<double> spacings) {
    double inverseSquaredSum = 0.0;
    double smallest = std::numeric_limits<double>::infinity();
    for (double spacing : spacings) {
        if (spacing <= 0.0) {
            continue;
        }
        inverseSquaredSum += 1.0 / (spacing * spacing);
        smallest = std::min(smallest, spacing);
    }

    double limit = std::numeric_limits<double>::infinity();

    // von Neumann: the forward-Euler amplification factor for the diffusive
    // term is 1 - 4 dt nu sum(1/h^2) in the worst mode, so staying inside
    // |g| <= 1 needs dt <= 1/(2 nu sum(1/h^2)). The factor of two is the one
    // TransientDiffusionSolver was missing.
    if (effectiveViscosity > 0.0 && inverseSquaredSum > 0.0) {
        limit = std::min(limit, 1.0 / (2.0 * effectiveViscosity * inverseSquaredSum));
    }
    // CFL: information may not cross more than one cell per step.
    if (velocityScale > 0.0 && smallest < std::numeric_limits<double>::infinity()) {
        limit = std::min(limit, smallest / velocityScale);
    }

    return kExplicitStabilitySafety * limit;
}

} // namespace aether::solver
