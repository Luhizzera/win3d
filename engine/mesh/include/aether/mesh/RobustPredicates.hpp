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
