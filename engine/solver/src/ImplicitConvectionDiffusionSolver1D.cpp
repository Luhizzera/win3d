#include "aether/solver/ImplicitConvectionDiffusionSolver1D.hpp"

#include "aether/solver/ConvectionLimiter.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace aether::solver {

ImplicitConvectionDiffusionSolver1D::ImplicitConvectionDiffusionSolver1D(std::size_t nx, double length,
                                                                           double velocity, double diffusivity,
                                                                           double source, double leftValue,
                                                                           double rightValue)
    : nx_(nx), length_(length), h_(length / static_cast<double>(nx)), velocity_(nx, velocity),
      diffusivity_(nx, diffusivity), source_(source), leftValue_(leftValue), rightValue_(rightValue),
      phi_(nx, 0.0) {}

void ImplicitConvectionDiffusionSolver1D::setVelocityField(const std::vector<double>& velocityPerCell) {
    if (velocityPerCell.size() != nx_) {
        throw std::invalid_argument("velocity field size must equal the cell count");
    }
    velocity_ = velocityPerCell;
    // Invalidates BOTH the cached matrix extraction and anything built
    // from it -- forgetting the second one would leave a stale ILU(0)
    // factorization of the previous coefficients in place.
    diag_.clear();
    iluDiag_.clear();
}

void ImplicitConvectionDiffusionSolver1D::setDiffusivityField(const std::vector<double>& diffusivityPerCell) {
    if (diffusivityPerCell.size() != nx_) {
        throw std::invalid_argument("diffusivity field size must equal the cell count");
    }
    diffusivity_ = diffusivityPerCell;
    diag_.clear();
    iluDiag_.clear();
}

void ImplicitConvectionDiffusionSolver1D::setPreconditioner(Preconditioner preconditioner) {
    preconditioner_ = preconditioner;
    iluDiag_.clear();
}

double ImplicitConvectionDiffusionSolver1D::velocityAtFace(std::size_t face) const {
    if (face == 0) {
        return velocity_[0];
    }
    if (face >= nx_) {
        return velocity_[nx_ - 1];
    }
    return 0.5 * (velocity_[face - 1] + velocity_[face]);
}

double ImplicitConvectionDiffusionSolver1D::diffusivityAtFace(std::size_t face) const {
    if (face == 0) {
        return diffusivity_[0];
    }
    if (face >= nx_) {
        return diffusivity_[nx_ - 1];
    }
    return 0.5 * (diffusivity_[face - 1] + diffusivity_[face]);
}

double ImplicitConvectionDiffusionSolver1D::cellCenterX(std::size_t i) const {
    return (static_cast<double>(i) + 0.5) * h_;
}

double ImplicitConvectionDiffusionSolver1D::phiAt(const std::vector<double>& x, long long i) const {
    if (i < 0) {
        return 2.0 * leftValue_ - x[0];
    }
    if (i >= static_cast<long long>(nx_)) {
        return 2.0 * rightValue_ - x[nx_ - 1];
    }
    return x[static_cast<std::size_t>(i)];
}

double ImplicitConvectionDiffusionSolver1D::stencilAt(const std::vector<double>& x, std::size_t i) const {
    const auto li = static_cast<long long>(i);
    const double phiC = x[i];
    const double phiE = phiAt(x, li + 1);
    const double phiW = phiAt(x, li - 1);

    // Face f sits between cells f-1 and f, so cell i's west face is i and
    // its east face is i+1.
    const double uE = velocityAtFace(i + 1);
    const double uW = velocityAtFace(i);
    const double gammaE = diffusivityAtFace(i + 1);
    const double gammaW = diffusivityAtFace(i);

    const double fluxE = uE >= 0.0 ? uE * phiC : uE * phiE;
    const double fluxW = uW >= 0.0 ? uW * phiW : uW * phiC;
    const double diffE = gammaE * (phiE - phiC) / h_;
    const double diffW = gammaW * (phiC - phiW) / h_;

    return (fluxE - fluxW) - (diffE - diffW);
}

std::vector<double> ImplicitConvectionDiffusionSolver1D::applyOperator(const std::vector<double>& x) const {
    const std::vector<double> zero(nx_, 0.0);
    std::vector<double> result(nx_);
    for (std::size_t i = 0; i < nx_; ++i) {
        result[i] = stencilAt(x, i) - stencilAt(zero, i);
    }
    return result;
}

std::vector<double> ImplicitConvectionDiffusionSolver1D::rightHandSide() const {
    const std::vector<double> zero(nx_, 0.0);
    std::vector<double> rhs(nx_);
    for (std::size_t i = 0; i < nx_; ++i) {
        rhs[i] = source_ * h_ - stencilAt(zero, i);
    }
    // Deferred correction: the part of the convective flux the upwind matrix
    // cannot represent, evaluated from the previous outer iterate. Empty
    // under FirstOrderUpwind, which is what keeps that path unchanged.
    if (!convectionCorrection_.empty()) {
        for (std::size_t i = 0; i < nx_; ++i) {
            rhs[i] += convectionCorrection_[i];
        }
    }
    return rhs;
}

double ImplicitConvectionDiffusionSolver1D::convectedFaceValue(const std::vector<double>& x,
                                                                std::size_t face) const {
    const auto li = static_cast<long long>(face);
    // Face `face` sits between cells face-1 and face.
    const double faceVelocity = velocityAtFace(face);
    const double phiLeft = phiAt(x, li - 1);
    const double phiRight = phiAt(x, li);
    const bool leftIsUpwind = faceVelocity >= 0.0;
    const double upwind = leftIsUpwind ? phiLeft : phiRight;
    const double downwind = leftIsUpwind ? phiRight : phiLeft;
    const double central = 0.5 * (phiLeft + phiRight);

    if (scheme_ == ConvectionScheme::FirstOrderUpwind) {
        return upwind;
    }
    if (scheme_ == ConvectionScheme::Central) {
        return central;
    }

    const double difference = downwind - upwind;
    if (faceDifferenceIsNegligible(difference, upwind, downwind)) {
        return upwind;
    }
    // The upwind-upwind cell, which a structured grid has and a tetrahedral
    // mesh does not -- see ConvectionLimiter.hpp for why the ratio built from
    // it is the same quantity the unstructured side builds from a gradient.
    const double farUpwind = leftIsUpwind ? phiAt(x, li - 2) : phiAt(x, li + 1);
    const double ratio = (upwind - farUpwind) / difference;
    return upwind + vanLeerLimiter(ratio) * (central - upwind);
}

std::vector<double> ImplicitConvectionDiffusionSolver1D::convectionCorrection(
    const std::vector<double>& x) const {
    std::vector<double> correction(nx_, 0.0);
    if (scheme_ == ConvectionScheme::FirstOrderUpwind) {
        return correction;
    }
    // Cell i's west face is i and its east face is i+1; the correction enters
    // with the opposite sign to the flux, since it moves to the other side of
    // the equation.
    for (std::size_t i = 0; i < nx_; ++i) {
        const double uE = velocityAtFace(i + 1);
        const double uW = velocityAtFace(i);
        const double phiE = convectedFaceValue(x, i + 1);
        const double phiW = convectedFaceValue(x, i);
        const double upwindE = uE >= 0.0 ? phiAt(x, static_cast<long long>(i))
                                          : phiAt(x, static_cast<long long>(i) + 1);
        const double upwindW = uW >= 0.0 ? phiAt(x, static_cast<long long>(i) - 1)
                                          : phiAt(x, static_cast<long long>(i));
        correction[i] = -((uE * phiE - uW * phiW) - (uE * upwindE - uW * upwindW));
    }
    return correction;
}

void ImplicitConvectionDiffusionSolver1D::outerSweeps(const std::function<void()>& innerSolve,
                                                       double tolerance) {
    if (scheme_ == ConvectionScheme::FirstOrderUpwind) {
        convectionCorrection_.clear();
        innerSolve();
        return;
    }
    constexpr std::size_t kMaxOuterSweeps = 300;
    for (std::size_t sweep = 0; sweep < kMaxOuterSweeps; ++sweep) {
        convectionCorrection_ = convectionCorrection(phi_);
        const std::vector<double> previous = phi_;
        innerSolve();
        double change = 0.0;
        for (std::size_t i = 0; i < nx_; ++i) {
            change = std::max(change, std::fabs(phi_[i] - previous[i]));
        }
        if (change < tolerance) {
            break;
        }
    }
}

double ImplicitConvectionDiffusionSolver1D::maxCellPeclet() const {
    double worst = 0.0;
    for (std::size_t face = 0; face <= nx_; ++face) {
        const double gamma = diffusivityAtFace(face);
        if (gamma > 0.0) {
            worst = std::max(worst, std::fabs(velocityAtFace(face)) * h_ / gamma);
        }
    }
    return worst;
}

double ImplicitConvectionDiffusionSolver1D::maxBoundednessViolation() const {
    const double low = std::min(leftValue_, rightValue_);
    const double high = std::max(leftValue_, rightValue_);
    const double range = high - low;
    if (range <= 0.0) {
        return 0.0;
    }
    double worst = 0.0;
    for (double value : phi_) {
        worst = std::max(worst, std::max(low - value, value - high) / range);
    }
    return std::max(worst, 0.0);
}

double ImplicitConvectionDiffusionSolver1D::dot(const std::vector<double>& a, const std::vector<double>& b) {
    double result = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

void ImplicitConvectionDiffusionSolver1D::extractTridiagonal() {
    if (!diag_.empty()) {
        return; // already valid; invalidated by any coefficient change
    }
    sub_.assign(nx_, 0.0);
    diag_.assign(nx_, 0.0);
    super_.assign(nx_, 0.0);

    std::vector<double> probe(nx_, 0.0);
    for (std::size_t j = 0; j < nx_; ++j) {
        probe[j] = 1.0;
        const std::vector<double> column = applyOperator(probe);
        probe[j] = 0.0;

        diag_[j] = column[j];
        if (j + 1 < nx_) {
            sub_[j + 1] = column[j + 1]; // row j+1's coefficient on phi[j]
        }
        if (j > 0) {
            super_[j - 1] = column[j - 1]; // row j-1's coefficient on phi[j]
        }

        // The ILU(0)-is-exact argument below depends on this operator
        // genuinely being tridiagonal, so check it rather than assume it.
        for (std::size_t row = 0; row < nx_; ++row) {
            const bool onBand = row + 1 == j || row == j || row == j + 1;
            if (!onBand && column[row] != 0.0) {
                throw std::logic_error("convection-diffusion operator is not tridiagonal");
            }
        }
    }
}

void ImplicitConvectionDiffusionSolver1D::buildIncompleteLU() {
    if (!iluDiag_.empty()) {
        return;
    }
    // ILU(0) on a tridiagonal matrix drops nothing (LU of a tridiagonal
    // matrix has no fill-in outside the existing band), so this is the
    // exact LU factorization -- the classic Thomas-algorithm recurrence.
    iluDiag_.assign(nx_, 0.0);
    iluSub_.assign(nx_, 0.0);
    iluDiag_[0] = diag_[0];
    for (std::size_t i = 1; i < nx_; ++i) {
        iluSub_[i] = sub_[i] / iluDiag_[i - 1];
        iluDiag_[i] = diag_[i] - iluSub_[i] * super_[i - 1];
    }
}

std::vector<double> ImplicitConvectionDiffusionSolver1D::applyPreconditioner(
    const std::vector<double>& v) const {
    switch (preconditioner_) {
    case Preconditioner::None:
        return v;
    case Preconditioner::Jacobi: {
        std::vector<double> result(nx_);
        for (std::size_t i = 0; i < nx_; ++i) {
            result[i] = v[i] / diag_[i];
        }
        return result;
    }
    case Preconditioner::IncompleteLU: {
        std::vector<double> y(nx_);
        y[0] = v[0];
        for (std::size_t i = 1; i < nx_; ++i) { // forward solve L y = v
            y[i] = v[i] - iluSub_[i] * y[i - 1];
        }
        std::vector<double> result(nx_);
        result[nx_ - 1] = y[nx_ - 1] / iluDiag_[nx_ - 1];
        for (std::size_t i = nx_ - 1; i-- > 0;) { // back solve U x = y
            result[i] = (y[i] - super_[i] * result[i + 1]) / iluDiag_[i];
        }
        return result;
    }
    }
    return v;
}

void ImplicitConvectionDiffusionSolver1D::prepareSolve() {
    extractTridiagonal();
    if (preconditioner_ == Preconditioner::IncompleteLU) {
        buildIncompleteLU();
    }
}

std::size_t ImplicitConvectionDiffusionSolver1D::solveGaussSeidel(std::size_t maxIterations,
                                                                   double tolerance) {
    std::size_t iterations = 0;
    outerSweeps([&] { iterations = solveGaussSeidelOnce(maxIterations, tolerance); }, tolerance);
    return iterations;
}

std::size_t ImplicitConvectionDiffusionSolver1D::solveGaussSeidelOnce(std::size_t maxIterations,
                                                                      double tolerance) {
    // stencilAt(x, i) is affine in x[i] alone (holding the rest fixed),
    // with slope exactly diag_[i] -- so the update that zeroes cell i's
    // residual, given every other cell's current value, is a plain
    // residual/diagonal correction. Using the *extracted* diagonal rather
    // than a hand-derived formula covers both upwind branches and the two
    // boundary rows (where the ghost mirror folds a neighbor coefficient
    // onto the diagonal) without special-casing any of them.
    extractTridiagonal();
    std::size_t iteration = 0;
    for (; iteration < maxIterations; ++iteration) {
        double maxChange = 0.0;
        for (std::size_t i = 0; i < nx_; ++i) {
            const double residual = stencilAt(phi_, i) - source_ * h_;
            const double delta = -residual / diag_[i];
            phi_[i] += delta;
            maxChange = std::max(maxChange, std::fabs(delta));
        }
        if (maxChange < tolerance) {
            ++iteration;
            break;
        }
    }
    return iteration;
}

long long ImplicitConvectionDiffusionSolver1D::solveBiCGStab(std::size_t maxIterations,
                                                              double tolerance) {
    long long iterations = 0;
    outerSweeps([&] { iterations = solveBiCGStabOnce(maxIterations, tolerance); }, tolerance);
    return iterations;
}

long long ImplicitConvectionDiffusionSolver1D::solveBiCGStabOnce(std::size_t maxIterations,
                                                                  double tolerance) {
    residualHistory_.clear();
    prepareSolve();
    // Left preconditioning: solve M^-1 A x = M^-1 b. With
    // Preconditioner::None both wrappers are the identity, so this reduces
    // exactly to the unpreconditioned algorithm.
    const auto operatorApply = [this](const std::vector<double>& x) {
        return applyPreconditioner(applyOperator(x));
    };
    const std::vector<double> b = applyPreconditioner(rightHandSide());
    const double bNorm = std::max(std::sqrt(dot(b, b)), 1e-300);

    std::vector<double> r(nx_);
    {
        const std::vector<double> ax = operatorApply(phi_);
        for (std::size_t i = 0; i < nx_; ++i) {
            r[i] = b[i] - ax[i];
        }
    }
    const std::vector<double> rHat = r; // arbitrary "shadow" residual, standard BiCGSTAB choice r0_hat = r0

    double rho = 1.0;
    double alpha = 1.0;
    double omega = 1.0;
    std::vector<double> v(nx_, 0.0);
    std::vector<double> p(nx_, 0.0);

    // **Breakdown detection must be relative, and a breakdown reached with
    // an already-small residual is convergence, not failure.** Both of
    // these were gotten wrong in this function's first version and were
    // caught by a test, not by review:
    //
    // (1) The original check compared these inner products against an
    //     absolute 1e-300. But rho and rHat.v are inner products of
    //     residual-scale vectors, so an absolute floor means something
    //     completely different from one problem to the next. The tests
    //     below compare against machine epsilon times the norms actually
    //     involved, which detects the real failure mode -- the two vectors
    //     having lost any usable component along each other -- rather than
    //     mere smallness. **The factor is epsilon SQUARED, not epsilon,
    //     and that distinction is not cosmetic** -- an intermediate version
    //     of this function used plain epsilon and it fired *spuriously*:
    //     rho legitimately decays to epsilon-scale as BiCGSTAB approaches
    //     convergence without anything having gone wrong, so an
    //     epsilon-scale test aborts healthy solves. Caught because the C++
    //     test reported a breakdown on a case the Python bindings solved in
    //     107 iterations -- identical inputs (verified bit-for-bit) and the
    //     identical library, differing only by which binary's
    //     floating-point path got there first, which is exactly the
    //     signature of a threshold sitting on a knife edge. The real
    //     failure being guarded against is division blowing up, which needs
    //     a near-exact zero, not a small number.
    //
    // (2) BiCGSTAB is well known to stagnate or break down as it
    //     approaches the accuracy a problem's conditioning allows.
    //     Reporting failure there is misleading when the iterate in hand is
    //     perfectly good. Measured on this very class: a
    //     variable-coefficient case asked for 1e-12, broke down on the
    //     rHat.v test, and was holding a relative residual of 1.69e-12 --
    //     an answer the first version threw away as a failure. Now every
    //     breakdown path re-measures the true residual and reports success
    //     if it meets the tolerance.
    const double epsilonSquared = std::numeric_limits<double>::epsilon() *
                                   std::numeric_limits<double>::epsilon();
    const double rHatNorm = std::sqrt(dot(rHat, rHat));

    const auto concludeBreakdown = [&](std::size_t iterationsDone) -> long long {
        const std::vector<double> ax = operatorApply(phi_);
        double residualNormSquared = 0.0;
        for (std::size_t i = 0; i < nx_; ++i) {
            const double value = b[i] - ax[i];
            residualNormSquared += value * value;
        }
        const double relativeResidual = std::sqrt(residualNormSquared) / bNorm;
        residualHistory_.push_back(relativeResidual);
        if (relativeResidual < tolerance) {
            return static_cast<long long>(iterationsDone);
        }
        return -1;
    };

    for (std::size_t iteration = 0; iteration < maxIterations; ++iteration) {
        const double rhoNew = dot(rHat, r);
        if (std::fabs(rhoNew) <= epsilonSquared * rHatNorm * std::sqrt(dot(r, r))) {
            return concludeBreakdown(iteration + 1);
        }
        const double beta = (rhoNew / rho) * (alpha / omega);
        for (std::size_t i = 0; i < nx_; ++i) {
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        v = operatorApply(p);
        const double rHatDotV = dot(rHat, v);
        if (std::fabs(rHatDotV) <= epsilonSquared * rHatNorm * std::sqrt(dot(v, v))) {
            return concludeBreakdown(iteration + 1);
        }
        alpha = rhoNew / rHatDotV;

        std::vector<double> s(nx_);
        for (std::size_t i = 0; i < nx_; ++i) {
            s[i] = r[i] - alpha * v[i];
        }
        if (std::sqrt(dot(s, s)) / bNorm < tolerance) {
            for (std::size_t i = 0; i < nx_; ++i) {
                phi_[i] += alpha * p[i];
            }
            residualHistory_.push_back(std::sqrt(dot(s, s)) / bNorm);
            return static_cast<long long>(iteration + 1);
        }

        const std::vector<double> t = operatorApply(s);
        const double tDotT = dot(t, t);
        if (tDotT <= epsilonSquared * epsilonSquared * dot(s, s)) {
            // A annihilating s means omega has no meaningful minimizer.
            return concludeBreakdown(iteration + 1);
        }
        omega = dot(t, s) / tDotT;

        for (std::size_t i = 0; i < nx_; ++i) {
            phi_[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }
        const double relativeResidual = std::sqrt(dot(r, r)) / bNorm;
        residualHistory_.push_back(relativeResidual);
        if (relativeResidual < tolerance) {
            return static_cast<long long>(iteration + 1);
        }
        rho = rhoNew;
    }
    return static_cast<long long>(maxIterations);
}

long long ImplicitConvectionDiffusionSolver1D::solveGmres(std::size_t restart,
                                                           std::size_t maxIterations,
                                                           double tolerance) {
    long long iterations = 0;
    outerSweeps([&] { iterations = solveGmresOnce(restart, maxIterations, tolerance); }, tolerance);
    return iterations;
}

long long ImplicitConvectionDiffusionSolver1D::solveGmresOnce(std::size_t restart, std::size_t maxIterations,
                                                            double tolerance) {
    residualHistory_.clear();
    prepareSolve();
    // Left preconditioning, same as solveBiCGStab(): with
    // Preconditioner::None both wrappers are the identity and this reduces
    // exactly to unpreconditioned GMRES.
    const auto operatorApply = [this](const std::vector<double>& x) {
        return applyPreconditioner(applyOperator(x));
    };
    const std::vector<double> b = applyPreconditioner(rightHandSide());
    const double bNorm = std::max(std::sqrt(dot(b, b)), 1e-300);
    const std::size_t m = std::max<std::size_t>(restart, 1);

    // Same relative-threshold reasoning as solveBiCGStab() above: these are
    // norms of Krylov vectors, so the meaningful comparison is against
    // machine epsilon times the scale they actually live at (the current
    // cycle's initial residual norm), not an absolute constant.
    const double epsilon = std::numeric_limits<double>::epsilon();
    std::size_t totalIterations = 0;

    while (totalIterations < maxIterations) {
        // --- Each restart cycle begins afresh from the current iterate,
        // which is why restarted GMRES stays monotone globally: the zero
        // correction is always available inside the new Krylov subspace,
        // so a cycle can never make the residual worse than it started.
        std::vector<double> r(nx_);
        {
            const std::vector<double> ax = operatorApply(phi_);
            for (std::size_t i = 0; i < nx_; ++i) {
                r[i] = b[i] - ax[i];
            }
        }
        const double beta = std::sqrt(dot(r, r));
        if (beta / bNorm < tolerance) {
            return static_cast<long long>(totalIterations);
        }

        std::vector<std::vector<double>> basis; // Arnoldi basis vectors v_0..v_k
        basis.reserve(m + 1);
        basis.push_back(r);
        for (double& value : basis[0]) {
            value /= beta;
        }

        // hessenberg[j] is column j of the (rotated) Hessenberg matrix,
        // holding its j+2 potentially-nonzero entries.
        std::vector<std::vector<double>> hessenberg;
        std::vector<double> cosines(m, 0.0);
        std::vector<double> sines(m, 0.0);
        std::vector<double> g(m + 1, 0.0); // rotated RHS of the least-squares problem
        g[0] = beta;

        std::size_t completed = 0;
        bool converged = false;
        bool breakdown = false;

        for (std::size_t j = 0; j < m && totalIterations < maxIterations; ++j) {
            std::vector<double> w = operatorApply(basis[j]);

            // Modified Gram-Schmidt: each projection is subtracted from the
            // *already-updated* w rather than the original. Mathematically
            // identical to classical Gram-Schmidt, numerically far better --
            // the standard choice inside Arnoldi for exactly this reason.
            std::vector<double> column(j + 2, 0.0);
            for (std::size_t i = 0; i <= j; ++i) {
                column[i] = dot(w, basis[i]);
                for (std::size_t n = 0; n < nx_; ++n) {
                    w[n] -= column[i] * basis[i][n];
                }
            }
            const double arnoldiNorm = std::sqrt(dot(w, w));
            column[j + 1] = arnoldiNorm;

            // Replay every previous rotation onto the new column, then
            // build the one that annihilates its subdiagonal entry -- this
            // keeps an up-to-date QR factorization of the Hessenberg matrix
            // incrementally, so |g[j+1]| below *is* the exact residual norm
            // of the least-squares solution, available without solving it.
            for (std::size_t i = 0; i < j; ++i) {
                const double rotated = cosines[i] * column[i] + sines[i] * column[i + 1];
                column[i + 1] = -sines[i] * column[i] + cosines[i] * column[i + 1];
                column[i] = rotated;
            }
            const double denominator = std::hypot(column[j], column[j + 1]);
            if (denominator <= epsilon * beta) {
                breakdown = true;
                break;
            }
            cosines[j] = column[j] / denominator;
            sines[j] = column[j + 1] / denominator;
            column[j] = denominator;
            column[j + 1] = 0.0;
            g[j + 1] = -sines[j] * g[j];
            g[j] = cosines[j] * g[j];

            hessenberg.push_back(std::move(column));
            ++completed;
            ++totalIterations;

            const double relativeResidual = std::fabs(g[j + 1]) / bNorm;
            residualHistory_.push_back(relativeResidual);

            // A zero Arnoldi norm is a "happy breakdown": the Krylov
            // subspace is already invariant, so the least-squares solution
            // over it is the exact solution -- success, not failure.
            if (relativeResidual < tolerance || arnoldiNorm <= epsilon * beta) {
                converged = true;
                break;
            }
            for (double& value : w) {
                value /= arnoldiNorm;
            }
            basis.push_back(std::move(w));
        }

        // Back-substitute the (upper-triangular, post-rotation) system
        // H*y = g, then apply the correction x += V*y. Done at the end of
        // every cycle, converged or not, so a restart resumes from the best
        // iterate this cycle found rather than discarding its work.
        if (completed > 0) {
            std::vector<double> y(completed, 0.0);
            for (std::size_t row = completed; row-- > 0;) {
                double sum = g[row];
                for (std::size_t col = row + 1; col < completed; ++col) {
                    sum -= hessenberg[col][row] * y[col];
                }
                if (std::fabs(hessenberg[row][row]) <= epsilon * beta) {
                    return -1;
                }
                y[row] = sum / hessenberg[row][row];
            }
            for (std::size_t col = 0; col < completed; ++col) {
                for (std::size_t n = 0; n < nx_; ++n) {
                    phi_[n] += y[col] * basis[col][n];
                }
            }
        }

        if (breakdown) {
            return -1;
        }
        if (converged) {
            return static_cast<long long>(totalIterations);
        }
    }
    return static_cast<long long>(totalIterations);
}

} // namespace aether::solver
