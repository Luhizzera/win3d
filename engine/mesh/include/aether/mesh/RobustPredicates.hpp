#pragma once

#include "aether/core/Vector3.hpp"

namespace aether::mesh {

// Exact-arithmetic replacements for the two geometric sign tests
// DelaunayTetrahedralization3D's Bowyer-Watson implementation depends on
// most heavily: orientation (which side of a plane a point is on) and
// in-sphere (Delaunay's own defining circumsphere-emptiness test).
//
// **Why this exists**: root-caused a real, if rare (~1/50 random small
// point sets), bug where tetrahedralize() produced a tetrahedralization
// with a genuine gap. Traced it to inCircumsphere()'s plain double-
// precision 4x4 determinant giving an unreliable answer for a near-
// degenerate ("sliver") input tetrahedron -- one specific failing case
// had a tetrahedron with only ~2.7% of the volume its own edge lengths
// would predict for a well-shaped tetrahedron, giving it a circumsphere
// radius wildly out of scale with the point cloud it belonged to. This is
// the standard, well-documented reason production Delaunay implementations
// (CGAL, TetGen, Qhull) use exact or adaptive-precision arithmetic for
// these specific predicates instead of naive floating-point determinants.
//
// **How this achieves exactness**: every IEEE754 double is *exactly*
// representable as significand * 2^exponent (a dyadic rational) -- no
// approximation is needed to convert one, only to compute with it in
// plain floating point afterward. Every input coordinate is converted to
// such an exact (aether::core::BigInt significand, int exponent) pair up
// front, and every subsequent operation (subtraction, multiplication,
// addition) is carried out on that exact representation -- add/subtract
// align to a common exponent via BigInt::shiftLeft() first, multiply adds
// exponents and multiplies significands -- so no floating-point rounding
// ever re-enters the computation. Only the final *sign* is extracted, not
// a numeric value, which is all either predicate actually needs.
// aether::core::BigInt's own header explains why its fixed (1280-bit)
// width is enough: these are degree-3 (orientation) and degree-5
// (in-sphere) polynomials in the input coordinates, needing at most
// roughly 270 bits for realistic double inputs, well under budget.
//
// Both formulas exactly mirror the plain-double versions that predate
// this (still present, e.g. in makeTetrahedron()'s own orientation check)
// -- same cofactor expansions, same sign conventions -- just routed
// through exact arithmetic instead of double arithmetic throughout.

// **A floating-point filter runs first, and it does not change any
// answer.** Exact BigInt arithmetic is correct but expensive, and
// DelaunayTetrahedralization3D calls inSphere3D once per tetrahedron per
// inserted point -- measured as the dominant cost of tetrahedralize()
// (see below). The filter computes the same determinant in plain doubles
// *together with a rigorous bound on its own accumulated roundoff*; when
// |determinant| exceeds that bound the double sign is provably correct
// and is returned directly, and only the remaining near-degenerate cases
// pay for exact arithmetic. This is the standard technique (Shewchuk's
// adaptive predicates; CGAL/TetGen/Qhull all do it) and is what the
// "exact or adaptive-precision" note above was always pointing at.
//
// The bound constants are deliberately several times more conservative
// than the tightest published ones: over-estimating the error only sends
// more cases down the exact path (slower, never wrong), while
// under-estimating it would silently return a wrong sign -- an asymmetry
// worth paying for. `orientation3DExact`/`inSphere3DExact` below expose
// the unfiltered path so the two can be compared directly, which is how
// the filter is validated (mesh_tests runs both over randomized and
// deliberately degenerate configurations and requires identical signs).

// Sign of (b-a) . ((c-a) x (d-a)): positive means (a,b,c,d) is a
// positively-oriented tetrahedron (matches every other orientation check
// in this project). +1 / 0 / -1.
int orientation3D(const core::Vector3& a, const core::Vector3& b, const core::Vector3& c,
                   const core::Vector3& d);

// Sign of the lifted-paraboloid 4x4 determinant testing whether p lies
// inside the circumsphere of tetrahedron abcd -- +1 means strictly
// inside, matching DelaunayTetrahedralization3D::inCircumsphere's
// existing convention (which assumes abcd is positively oriented, same
// requirement carried over here). +1 / 0 / -1.
int inSphere3D(const core::Vector3& a, const core::Vector3& b, const core::Vector3& c,
                const core::Vector3& d, const core::Vector3& p);

// The unfiltered, always-exact paths the two functions above fall back
// to. Identical results by construction -- the filter only ever decides
// whether the exact computation can be *skipped*, never what the answer
// is -- so these exist for validating that claim rather than for callers
// to choose between. Prefer the filtered versions everywhere else: they
// are the same answer for less work.
int orientation3DExact(const core::Vector3& a, const core::Vector3& b, const core::Vector3& c,
                        const core::Vector3& d);
int inSphere3DExact(const core::Vector3& a, const core::Vector3& b, const core::Vector3& c,
                     const core::Vector3& d, const core::Vector3& p);

// The 2D counterparts, for DelaunayTriangulation2D and
// PolygonTriangulation2D -- same exact-arithmetic approach, one dimension
// down. Both ignore z (every point in those classes keeps z = 0 by
// construction).

// Sign of the z-component of (b-a) x (c-a): positive means (a,b,c) winds
// counter-clockwise (matches the 2D classes' own CCW convention). Also
// exactly the sign of twice the signed area, so 0 means collinear.
// +1 / 0 / -1.
int orientation2D(const core::Vector3& a, const core::Vector3& b, const core::Vector3& c);

// Sign of the lifted-parabola 3x3 determinant testing whether p lies
// inside the circumcircle of triangle abc -- +1 means strictly inside,
// matching DelaunayTriangulation2D::inCircumcircle's existing convention
// (which assumes abc is CCW). +1 / 0 / -1.
int inCircle2D(const core::Vector3& a, const core::Vector3& b, const core::Vector3& c,
                const core::Vector3& p);

} // namespace aether::mesh
