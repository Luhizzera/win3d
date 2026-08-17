#pragma once

#include <cmath>

namespace aether::solver {

// The limiter itself, with no geometry attached, so the structured and
// unstructured solvers share one definition instead of two.
//
// **The shape of the blend.** A face value written as
//
//   phi_f = phi_C + psi * (phi_central - phi_C)
//
// is exactly upwind at psi = 0 and exactly the interpolated (second-order)
// value at psi = 1, so the limiter's whole job is choosing a point on that
// line. `phi_central` is whatever second-order interpolation the caller's
// geometry provides -- a plain average on a uniform grid, the
// distance-weighted one on a skewed mesh -- which is the only thing the two
// call sites need to differ about.
//
// **The ratio, and why the two callers' formulas are the same formula.**
// The classic structured ratio is r = (phi_C - phi_CC)/(phi_D - phi_C),
// comparing the upwind-side slope against the one being spanned; it needs an
// upwind-upwind cell, which a tetrahedral mesh does not have. The
// unstructured form r = 2 (grad(phi)_C . d)/(phi_D - phi_C) - 1 uses the
// upwind cell's own gradient instead. They are not two conventions: on a
// uniform grid the central-difference gradient is
// grad(phi)_C . d = (phi_D - phi_CC)/2, and substituting gives
//
//   2 * (phi_D - phi_CC)/2 / (phi_D - phi_C) - 1 = (phi_C - phi_CC)/(phi_D - phi_C)
//
// -- the classic ratio exactly. So both callers pass the same r and this
// function needs to know nothing about which grid it came from.
//
// van Leer's limiter, chosen smooth on purpose: r = 1 (a smooth field) gives
// psi = 1 and full second order, r <= 0 (an extremum) gives psi = 0 and
// upwind, and unlike minmod or superbee it has no kink -- a limiter with a
// kink makes the residual of a steady iteration non-differentiable and can
// stall it short of convergence.
inline double vanLeerLimiter(double ratio) {
    return (ratio + std::fabs(ratio)) / (1.0 + std::fabs(ratio));
}

// True when the difference across the face is too small for the ratio above
// to mean anything. A vanishing difference is not an extremum but a locally
// flat field, where every scheme agrees, so falling back to upwind there
// costs no accuracy. Scaled by the field's own magnitude so the test stays
// relative rather than depending on the units.
inline bool faceDifferenceIsNegligible(double difference, double upwindValue, double downwindValue) {
    const double scale = std::fabs(upwindValue) + std::fabs(downwindValue);
    return std::fabs(difference) <= 1e-12 * (scale + 1.0);
}

} // namespace aether::solver
