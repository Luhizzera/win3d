#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace aether::analysis {

// Module 12's "analytical" sub-area (first pass): automated, rule-based
// diagnostics computed over a solver's already-computed field data --
// deliberately not machine learning of any kind, just exact arithmetic
// derived from each quantity's own definition. Every number here is either
// a plain statistic or something derivable on paper (Courant number is a
// definition, not a recalled fact; the checkerboard index below is defined
// in this header, not borrowed from a table), so there is no literature-
// recall risk the way e.g. quoting a Ghia et al. benchmark value would be.

struct FieldStatistics {
    double minValue;
    double maxValue;
    double mean;
};

// Throws std::invalid_argument if `field` is empty.
FieldStatistics computeStatistics(const std::vector<double>& field);

// Multi-dimensional Courant (CFL) number at each cell,
// dt * (|u|/dx + |v|/dy), maximized over all cells -- the standard
// sum-form bound for explicit time-stepping stability under 2D advection
// (the same quantity `stableTimeStep()` in every 2D/3D Navier-Stokes
// solver in this project is implicitly designed to keep at or below 1;
// this function lets a caller measure it directly from a solved field
// instead of trusting that by construction). Throws std::invalid_argument
// if `u`/`v` sizes disagree.
double maxCourantNumber(const std::vector<double>& u, const std::vector<double>& v, double dx, double dy, double dt);

// A novel diagnostic (defined here, not borrowed from any external source):
// how strongly a 2D structured field exhibits the classic collocated-grid
// "checkerboard" (odd-even decoupling) pattern -- exactly the known
// weakness this project's own Navier-Stokes solvers document (no
// Rhie-Chow interpolation; see TaylorGreenVortexSolver2D/
// LidDrivenCavitySolver2D's class comments) but never actually measured
// until now.
//
// At each interior cell, comparing the cell's own value to the average of
// its 4 face-neighbours gives a residual that is exactly 0 for any field
// that is locally linear (in particular for any constant field) and, for a
// perfect single-cell checkerboard field(i,j) = A*(-1)^(i+j), is exactly
// 2*field(i,j) at every interior cell (each neighbour has the opposite
// sign, so their average is -field(i,j)). Normalizing the RMS of that
// residual by 2x the field's own overall RMS variation therefore gives
// **exactly 0** for a smooth/constant field and **exactly 1** for a
// perfect checkerboard, by construction -- both are checked directly in
// this module's own tests, not asserted without proof.
//
// Returns 0.0 for a field with zero variance (uniform field -- there is no
// meaningful "how checkerboarded" question to ask about a constant).
// Throws std::invalid_argument if nx < 3 or ny < 3 (no interior cells to
// measure) or if field.size() != nx*ny.
double checkerboardIndex(const std::vector<double>& field, std::size_t nx, std::size_t ny);

// One formatted line: "<name>: min=... max=... mean=... (n=...)".
std::string summarizeField(const std::string& name, const std::vector<double>& field);

} // namespace aether::analysis
