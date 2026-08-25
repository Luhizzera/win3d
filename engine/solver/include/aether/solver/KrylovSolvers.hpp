#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <vector>

namespace aether::solver {

// Matrix-free Krylov methods, in a header of their own.
//
// **Why they live here and not on a solver base.** GMRES was first written
// inside ImplicitConvectionDiffusionSolver1D and then written a second
// time as a member of UnstructuredFvmBase, which needed it for the
// pressure-velocity coupling of DIVIDA_TECNICA.md 4.3. That second copy
// recorded its own duplication honestly rather than hiding it -- and the
// fix is not simply to have one call the other, because neither is the
// right owner: a matrix-free Krylov method knows nothing about tetrahedral
// meshes or about one-dimensional convection-diffusion. It is linear
// algebra. Putting it in a base class of unstructured finite-volume
// solvers would have made every other caller inherit from something
// unrelated to reach it.
//
// This is the same lesson item 2.1 recorded when UnstructuredFvmBase was
// itself extracted, one level up: **porting a correction N times is the
// expensive version of the problem** -- and the place a shared thing lives
// should be decided by what it is, not by which caller happened to need it
// first.

namespace detail {

inline double krylovDot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

} // namespace detail

// **Restarted GMRES(m), matrix-free**, for operators Conjugate Gradient
// cannot solve because they are not symmetric positive definite.
//
// Modified Gram-Schmidt Arnoldi plus incremental Givens rotations: the
// rotations maintain an up-to-date QR factorization of the Hessenberg
// matrix, so the least-squares residual norm is available every iteration
// without solving the least-squares problem.
//
// `x` is both the initial guess and the result. Returns the number of
// matrix-vector products used.
//
// `onResidual`, if given, is called once per iteration with the relative
// residual -- the observability a caller needs to plot or assert on
// convergence history, without this function taking a position on what
// should be done with it.
//
// `brokeDown`, if given, is set true when the iteration stopped because a
// Givens denominator or a back-substitution pivot vanished. That is
// deliberately *reported* rather than folded into the return value:
// "converged", "ran out of iterations" and "broke down" are three
// different outcomes, and different callers reasonably care about
// different ones. A caller that only wants an answer can ignore it and
// check the residual itself, exactly as it would for CG.
template <typename ApplyOperator>
std::size_t gmres(ApplyOperator&& apply, const std::vector<double>& rhs, std::vector<double>& x,
                   std::size_t restart, std::size_t maxIterations, double tolerance,
                   const std::function<void(double)>& onResidual = {},
                   bool* brokeDown = nullptr) {
    using detail::krylovDot;

    const std::size_t n = rhs.size();
    const double bNorm = std::max(std::sqrt(krylovDot(rhs, rhs)), 1e-300);
    const std::size_t m = std::max<std::size_t>(restart, 1);
    const double epsilon = std::numeric_limits<double>::epsilon();
    std::size_t totalIterations = 0;
    if (brokeDown != nullptr) {
        *brokeDown = false;
    }

    while (totalIterations < maxIterations) {
        // --- Each restart cycle begins afresh from the current iterate,
        // which is why restarted GMRES stays monotone globally: the zero
        // correction is always available inside the new Krylov subspace,
        // so a cycle can never make the residual worse than it started.
        std::vector<double> r(n);
        {
            const std::vector<double> ax = apply(x);
            for (std::size_t i = 0; i < n; ++i) {
                r[i] = rhs[i] - ax[i];
            }
        }
        const double beta = std::sqrt(krylovDot(r, r));
        if (beta / bNorm < tolerance) {
            return totalIterations;
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
            std::vector<double> w = apply(basis[j]);

            // Modified Gram-Schmidt: each projection is subtracted from
            // the *already-updated* w rather than the original --
            // mathematically identical to classical Gram-Schmidt,
            // numerically far better, the standard choice inside Arnoldi
            // for exactly this reason.
            std::vector<double> column(j + 2, 0.0);
            for (std::size_t i = 0; i <= j; ++i) {
                column[i] = krylovDot(w, basis[i]);
                for (std::size_t k = 0; k < n; ++k) {
                    w[k] -= column[i] * basis[i][k];
                }
            }
            const double arnoldiNorm = std::sqrt(krylovDot(w, w));
            column[j + 1] = arnoldiNorm;

            // Replay every previous rotation onto the new column, then
            // build the one that annihilates its subdiagonal entry.
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
            if (onResidual) {
                onResidual(relativeResidual);
            }

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
            bool solvable = true;
            for (std::size_t row = completed; row-- > 0;) {
                double sum = g[row];
                for (std::size_t col = row + 1; col < completed; ++col) {
                    sum -= hessenberg[col][row] * y[col];
                }
                if (std::fabs(hessenberg[row][row]) <= epsilon * beta) {
                    solvable = false;
                    break;
                }
                y[row] = sum / hessenberg[row][row];
            }
            if (!solvable) {
                if (brokeDown != nullptr) {
                    *brokeDown = true;
                }
                return totalIterations; // keep the best x found so far
            }
            for (std::size_t col = 0; col < completed; ++col) {
                for (std::size_t k = 0; k < n; ++k) {
                    x[k] += y[col] * basis[col][k];
                }
            }
        }

        if (breakdown) {
            if (brokeDown != nullptr) {
                *brokeDown = true;
            }
            return totalIterations;
        }
        if (converged) {
            return totalIterations;
        }
    }
    return totalIterations;
}

} // namespace aether::solver
