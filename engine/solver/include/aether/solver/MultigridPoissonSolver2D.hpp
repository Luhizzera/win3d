#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 5, "next steps" list (Multigrid/GMRES/BiCGSTAB/preconditioners):
// geometric multigrid (V-cycle, correction scheme) for the 2D Poisson
// equation -nabla^2(phi) = source on a rectangle with Dirichlet boundary
// conditions -- the same equation SteadyDiffusionSolver::solve() and
// ::solveConjugateGradient() already solve, now with a solver whose
// iteration count is *grid-independent* (multigrid's classic selling
// point) rather than growing with resolution the way Gauss-Seidel and,
// to a lesser extent, Conjugate Gradient do.
//
// Deliberately scoped to a standalone 2D square-friendly solver rather
// than folded into DiffusionProblem/SteadyDiffusionSolver: coarsening
// requires nx and ny to each be a power of 2 (checked informally via
// documentation, not enforced at runtime -- passing a non-power-of-2 size
// will simply stop coarsening early, degrading to fewer V-cycle levels
// rather than crashing), which the general 6-face/anisotropic
// DiffusionProblem hierarchy was never designed around.
//
// Boundary treatment uses Dirichlet ghost-mirroring (ghost = 2*wallValue -
// interior), the same technique already validated in
// LidDrivenCavitySolver2D/MixingLengthLidDrivenCavitySolver2D -- chosen
// specifically *because* every cell in the solution array is then a
// genuine unknown (no boundary-fixed cell layer baked into the array, the
// convention SteadyDiffusionSolver uses instead), which makes coarsening
// cleanly pair 2x2 fine blocks into one coarse cell with no special-casing
// at the domain edges.
//
// V-cycle (correction scheme): pre-smooth with weighted Gauss-Seidel,
// compute the residual, restrict it to the coarse grid (simple average of
// each 2x2 fine block -- the natural cell-centered/finite-volume
// restriction for uniform cells) as the coarse grid's right-hand side for
// a *correction* equation (homogeneous, zero Dirichlet boundaries -- the
// fine grid's real boundary values are already exactly satisfied, only the
// interior error needs correcting), recurse, prolongate the coarse
// correction back (piecewise-constant injection: each fine cell in a
// block receives its coarse cell's correction -- simpler than bilinear
// interpolation, a documented simplification for this first cut), add it
// in, post-smooth. The coarsest level is smoothed many times as an
// approximate direct solve rather than solved exactly.
class MultigridPoissonSolver2D {
public:
    enum class Face { XMin, XMax, YMin, YMax };

    MultigridPoissonSolver2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY);

    void setBoundaryValue(Face face, double value);
    void setSourceTerm(double source);

    // Runs V-cycles until the finest level's residual L2 norm drops below
    // tolerance or maxVCycles is reached. Returns the number of V-cycles
    // actually run.
    std::size_t solve(std::size_t maxVCycles = 100, double tolerance = 1e-9, std::size_t preSweeps = 2,
                       std::size_t postSweeps = 2);

    double value(std::size_t i, std::size_t j) const;

private:
    struct Level {
        std::size_t nx;
        std::size_t ny;
        double hx;
        double hy;
        std::vector<double> phi;
        std::vector<double> rhs;
        double faceValue[4] = {0.0, 0.0, 0.0, 0.0}; // indexed by Face
    };

    std::size_t index(const Level& level, std::size_t i, std::size_t j) const { return i + j * level.nx; }
    double dirichletAt(const Level& level, std::size_t i, std::size_t j, int di, int dj) const;

    void smooth(Level& level, std::size_t sweeps) const;
    std::vector<double> computeResidual(const Level& level) const;
    static double residualNorm(const std::vector<double>& residual);
    std::vector<double> restrictToCoarse(const std::vector<double>& fineResidual, const Level& fineLevel,
                                          const Level& coarseLevel) const;
    void prolongateAndCorrect(Level& fineLevel, const Level& coarseLevel) const;

    void vCycle(std::size_t levelIndex, std::size_t preSweeps, std::size_t postSweeps);

    std::vector<Level> levels_;
};

} // namespace aether::solver
