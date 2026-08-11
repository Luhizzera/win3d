#include "aether/solver/MultigridPoissonSolver2D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

MultigridPoissonSolver2D::MultigridPoissonSolver2D(std::size_t nx, std::size_t ny, double lengthX,
                                                     double lengthY) {
    // Build the level hierarchy once, coarsening by 2x in each direction
    // until either dimension would drop below 4 cells.
    std::size_t curNx = nx;
    std::size_t curNy = ny;
    double curHx = lengthX / static_cast<double>(nx);
    double curHy = lengthY / static_cast<double>(ny);

    while (true) {
        Level level;
        level.nx = curNx;
        level.ny = curNy;
        level.hx = curHx;
        level.hy = curHy;
        level.phi.assign(curNx * curNy, 0.0);
        level.rhs.assign(curNx * curNy, 0.0);
        levels_.push_back(std::move(level));

        if (curNx % 2 != 0 || curNy % 2 != 0 || curNx / 2 < 4 || curNy / 2 < 4) {
            break;
        }
        curNx /= 2;
        curNy /= 2;
        curHx *= 2.0;
        curHy *= 2.0;
    }
}

void MultigridPoissonSolver2D::setBoundaryValue(Face face, double value) {
    levels_[0].faceValue[static_cast<int>(face)] = value;
}

void MultigridPoissonSolver2D::setSourceTerm(double source) {
    std::fill(levels_[0].rhs.begin(), levels_[0].rhs.end(), source);
}

double MultigridPoissonSolver2D::value(std::size_t i, std::size_t j) const {
    return levels_[0].phi[index(levels_[0], i, j)];
}

double MultigridPoissonSolver2D::dirichletAt(const Level& level, std::size_t i, std::size_t j, int di,
                                              int dj) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    if (ni < 0) {
        return 2.0 * level.faceValue[static_cast<int>(Face::XMin)] - level.phi[index(level, 0, j)];
    }
    if (ni >= static_cast<long long>(level.nx)) {
        return 2.0 * level.faceValue[static_cast<int>(Face::XMax)] -
               level.phi[index(level, level.nx - 1, j)];
    }
    if (nj < 0) {
        return 2.0 * level.faceValue[static_cast<int>(Face::YMin)] - level.phi[index(level, i, 0)];
    }
    if (nj >= static_cast<long long>(level.ny)) {
        return 2.0 * level.faceValue[static_cast<int>(Face::YMax)] -
               level.phi[index(level, i, level.ny - 1)];
    }
    return level.phi[index(level, static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))];
}

void MultigridPoissonSolver2D::smooth(Level& level, std::size_t sweeps) const {
    const double ax = 1.0 / (level.hx * level.hx);
    const double ay = 1.0 / (level.hy * level.hy);
    const double denom = 2.0 * ax + 2.0 * ay;
    for (std::size_t sweep = 0; sweep < sweeps; ++sweep) {
        for (std::size_t j = 0; j < level.ny; ++j) {
            for (std::size_t i = 0; i < level.nx; ++i) {
                const double left = dirichletAt(level, i, j, -1, 0);
                const double right = dirichletAt(level, i, j, 1, 0);
                const double down = dirichletAt(level, i, j, 0, -1);
                const double up = dirichletAt(level, i, j, 0, 1);
                const std::size_t idx = index(level, i, j);
                level.phi[idx] = (ax * (left + right) + ay * (down + up) + level.rhs[idx]) / denom;
            }
        }
    }
}

std::vector<double> MultigridPoissonSolver2D::computeResidual(const Level& level) const {
    const double ax = 1.0 / (level.hx * level.hx);
    const double ay = 1.0 / (level.hy * level.hy);
    const double denom = 2.0 * ax + 2.0 * ay;
    std::vector<double> residual(level.nx * level.ny);
    for (std::size_t j = 0; j < level.ny; ++j) {
        for (std::size_t i = 0; i < level.nx; ++i) {
            const double left = dirichletAt(level, i, j, -1, 0);
            const double right = dirichletAt(level, i, j, 1, 0);
            const double down = dirichletAt(level, i, j, 0, -1);
            const double up = dirichletAt(level, i, j, 0, 1);
            const std::size_t idx = index(level, i, j);
            const double appliedOperator = denom * level.phi[idx] - ax * (left + right) - ay * (down + up);
            residual[idx] = level.rhs[idx] - appliedOperator;
        }
    }
    return residual;
}

double MultigridPoissonSolver2D::residualNorm(const std::vector<double>& residual) {
    double sumSquares = 0.0;
    for (double r : residual) {
        sumSquares += r * r;
    }
    return std::sqrt(sumSquares / static_cast<double>(residual.size()));
}

std::vector<double> MultigridPoissonSolver2D::restrictToCoarse(const std::vector<double>& fineResidual,
                                                                 const Level& fineLevel,
                                                                 const Level& coarseLevel) const {
    std::vector<double> coarseRhs(coarseLevel.nx * coarseLevel.ny, 0.0);
    for (std::size_t cj = 0; cj < coarseLevel.ny; ++cj) {
        for (std::size_t ci = 0; ci < coarseLevel.nx; ++ci) {
            const std::size_t fi = 2 * ci;
            const std::size_t fj = 2 * cj;
            const double sum = fineResidual[index(fineLevel, fi, fj)] +
                                fineResidual[index(fineLevel, fi + 1, fj)] +
                                fineResidual[index(fineLevel, fi, fj + 1)] +
                                fineResidual[index(fineLevel, fi + 1, fj + 1)];
            coarseRhs[index(coarseLevel, ci, cj)] = sum / 4.0;
        }
    }
    return coarseRhs;
}

void MultigridPoissonSolver2D::prolongateAndCorrect(Level& fineLevel, const Level& coarseLevel) const {
    for (std::size_t cj = 0; cj < coarseLevel.ny; ++cj) {
        for (std::size_t ci = 0; ci < coarseLevel.nx; ++ci) {
            const double correction = coarseLevel.phi[index(coarseLevel, ci, cj)];
            const std::size_t fi = 2 * ci;
            const std::size_t fj = 2 * cj;
            fineLevel.phi[index(fineLevel, fi, fj)] += correction;
            fineLevel.phi[index(fineLevel, fi + 1, fj)] += correction;
            fineLevel.phi[index(fineLevel, fi, fj + 1)] += correction;
            fineLevel.phi[index(fineLevel, fi + 1, fj + 1)] += correction;
        }
    }
}

void MultigridPoissonSolver2D::vCycle(std::size_t levelIndex, std::size_t preSweeps,
                                       std::size_t postSweeps) {
    Level& level = levels_[levelIndex];

    if (levelIndex + 1 == levels_.size()) {
        smooth(level, preSweeps * 20); // coarsest level: approximate direct solve
        return;
    }

    smooth(level, preSweeps);
    const std::vector<double> residual = computeResidual(level);

    Level& coarse = levels_[levelIndex + 1];
    coarse.rhs = restrictToCoarse(residual, level, coarse);
    std::fill(coarse.phi.begin(), coarse.phi.end(), 0.0); // correction scheme: start each V-cycle at 0
    // Coarse levels always use homogeneous (zero) boundaries: the fine
    // grid's real Dirichlet values are already satisfied exactly, only the
    // interior error is being corrected here.

    vCycle(levelIndex + 1, preSweeps, postSweeps);

    prolongateAndCorrect(level, coarse);
    smooth(level, postSweeps);
}

std::size_t MultigridPoissonSolver2D::solve(std::size_t maxVCycles, double tolerance,
                                             std::size_t preSweeps, std::size_t postSweeps) {
    std::size_t cycle = 0;
    for (; cycle < maxVCycles; ++cycle) {
        vCycle(0, preSweeps, postSweeps);
        const double norm = residualNorm(computeResidual(levels_[0]));
        if (norm < tolerance) {
            ++cycle;
            break;
        }
    }
    return cycle;
}

} // namespace aether::solver
