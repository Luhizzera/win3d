#pragma once

#include "aether/core/Vector3.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace aether::mesh {

// Module 3: 3D Delaunay tetrahedralization of a point set, via the same
// Bowyer-Watson incremental algorithm already used by
// DelaunayTriangulation2D, generalized one dimension up: circumcircle ->
// circumsphere, super triangle -> super tetrahedron, and the star-shaped
// hole left by removing "bad" tetrahedra is now bounded by triangular
// faces (each belonging to exactly one bad tetrahedron) instead of edges.
//
// Points are genuine 3D points (unlike DelaunayTriangulation2D, which
// keeps z=0). Tetrahedra are always wound with positive orientation
// (signed volume > 0).
//
// This covers unconstrained tetrahedralization of a point cloud, plus
// Steiner point insertion (see insertSteinerPoint() below) for mesh
// refinement -- both interior points and points that extend the current
// convex hull. Boundary-constrained tetrahedralization (holes, multiple
// contours, preserving an imported surface mesh's own faces exactly) is a
// genuinely harder next step on this same core. recoverFacets() below
// covers two practical slices of it: the plain "usually already there"
// case, and (via tryFlipCoplanarQuadDiagonal()) the coplanar-quad diagonal
// ambiguity a convex boundary's own flat faces produce -- which is enough
// for any convex or mildly non-convex domain, a cylinder or box included.
// What is *not* covered, and remains genuinely open research, is a
// boundary with actual reflex (non-convex) structure severe enough that no
// local combinatorial swap or bounded Steiner subdivision recovers it at
// all -- the classic Schonhardt-polyhedron obstruction, where a
// tetrahedralization using only the boundary's own vertices provably does
// not exist for some inputs. No such input has ever actually been
// constructed or hit in this codebase; it is cited as the theoretical
// ceiling on what "guaranteed recovery" could mean, not an observed bug.
//
// A real bug once lived here and is now FIXED -- worth recording because
// the first diagnosis of it was wrong, and the wrong diagnosis was
// plausible. Symptom: tetrahedralize() rarely (~1/40 random small point
// sets) produced a tetrahedralization missing one thin, valid,
// hull-adjacent tetrahedron.
//
// The *incorrect* first diagnosis blamed floating-point precision in the
// in-sphere predicate, since the failing case did involve a near-degenerate
// "sliver" tetrahedron -- a genuinely suspicious-looking detail, and the
// standard reason production libraries (CGAL, TetGen, Qhull) use exact or
// adaptive-precision predicates. That led to building exact-arithmetic
// predicates (see RobustPredicates.hpp), which this class now uses.
//
// The *actual* cause, found by re-running the whole algorithm in exact
// rational arithmetic (Python Fractions -- no floating point anywhere) and
// seeing the identical gap persist: the super-tetrahedron was simply too
// small. See tetrahedralize()'s own comment for the full mechanism. The
// exact-arithmetic predicates were therefore not what fixed this bug --
// they are kept as an independently worthwhile robustness improvement,
// but the fix was enlarging the super-tetrahedron.
//
// Lesson worth keeping: "a sliver is involved, so it must be a precision
// problem" was a reasonable hypothesis that happened to be wrong. Re-doing
// the computation in exact arithmetic is what settled it -- a cheap,
// decisive test that should come *before* building infrastructure premised
// on the precision theory being right.
class DelaunayTetrahedralization3D {
public:
    struct Tetrahedron {
        std::array<std::size_t, 4> vertices;
    };

    // Adds a point (before tetrahedralize() is called) and returns its
    // index.
    std::size_t addPoint(double x, double y, double z);

    // Computes the Delaunay tetrahedralization of every point added so
    // far: insert points one at a time, remove every tetrahedron whose
    // circumsphere contains the new point, and re-tetrahedralize the
    // resulting cavity by connecting the new point to each boundary face.
    // A large enclosing "super tetrahedron" seeds the process and is
    // discarded at the end.
    void tetrahedralize();

    // Inserts one additional point into an already-computed
    // tetrahedralization (call after tetrahedralize()), updating the mesh
    // incrementally via the same Bowyer-Watson cavity-replacement step
    // tetrahedralize()'s own loop uses per point -- without rebuilding
    // from a super-tetrahedron. This is the standard way a real mesh
    // generator adds Steiner points for quality refinement (e.g. a
    // circumcenter, an edge midpoint, or any other point chosen to improve
    // tetrahedron shape) without retetrahedralizing from scratch.
    //
    // Supports points anywhere -- inside the current convex hull (the
    // ordinary Delaunay cavity-insertion case) or outside it (the hull
    // extends to include the new point, via the same "ghost tetrahedron"
    // visibility test findCavity() uses internally; see its own comment).
    // A degenerate case (e.g. p exactly coincident with an existing point,
    // or an exact coplanarity/co-sphericity tie) is still rejected
    // defensively -- every new tetrahedron the insertion would create must
    // be non-degenerate (nonzero volume) -- rather than silently producing
    // an invalid tetrahedralization.
    //
    // Returns true and appends the point (queryable via point(pointCount()-1))
    // if the insertion succeeded; returns false and leaves the mesh
    // unchanged otherwise.
    bool insertSteinerPoint(double x, double y, double z);

    // Boundary-constrained tetrahedralization, task #74's tractable slice.
    // Full generality (arbitrary non-convex boundary, guaranteed recovery
    // past the Schonhardt-polyhedron obstruction) is genuinely research-
    // level -- deliberately not attempted. What's here instead is the
    // practical subset production mesh generators (TetGen etc.) also lean
    // on for well-conditioned input: most boundary facets of a "nice"
    // domain (convex-ish, no adversarial near-degeneracies) are already
    // present in the *unconstrained* Delaunay tetrahedralization once every
    // boundary vertex has been added as an ordinary point; the few that
    // aren't can usually be recovered with a bounded number of Steiner
    // points. A facet given as three point indices is a boundary triangle
    // that must appear as some tetrahedron's face in the output.

    // Facets from `facets` NOT currently present as a face of any
    // tetrahedron (checked in either winding order). O(facets *
    // tetrahedronCount); a diagnostic/testing query, not a runtime one.
    std::vector<std::array<std::size_t, 3>> missingFacets(
        const std::vector<std::array<std::size_t, 3>>& facets) const;

    struct FacetRecoveryResult {
        // Elementary triangles now guaranteed present in the mesh: every
        // already-present input facet unchanged, plus -- for any facet
        // that needed help -- its recursive centroid-Steiner-point
        // subdivision. Feed this list to removeRegion() as the wall set.
        std::vector<std::array<std::size_t, 3>> recoveredFacets;
        // Input facets where even recursive centroid subdivision, down to
        // maxRounds levels deep, never produced a present sub-triangle.
        // Not expected for well-conditioned input, not eliminated in
        // principle (see the class header's Schonhardt-polyhedron note).
        std::vector<std::array<std::size_t, 3>> unrecovered;
    };

    // Attempts to recover every facet in `facets` that missingFacets()
    // would flag, in two phases -- combinatorial first, Steiner-point
    // subdivision second -- mirroring the order production mesh generators
    // (TetGen et al.) use:
    //
    //  1. tryFlipCoplanarQuadDiagonal(): for each pair of still-missing
    //     facets that together describe one half-quad each of the same
    //     coplanar quadrilateral, swaps which diagonal the existing
    //     tetrahedra use -- no new point, just two tetrahedra replaced by
    //     two others covering the identical volume. This is the technique
    //     DIVIDA_TECNICA.md/ROADMAP.md name as still missing (the "caso do
    //     quadrilátero coplanar"), and it is what fixes the one concretely
    //     reproduced failure this class has ever had: an axis-aligned
    //     cube's square faces, each an arbitrary choice between its two
    //     diagonals.
    //  2. Whatever the first phase couldn't reach: the original heuristic,
    //     inserting a missing facet's centroid as a Steiner point
    //     (insertSteinerPoint()) and recursing on the resulting three
    //     sub-triangles, up to maxRounds levels deep.
    //
    // Neither phase is a guarantee -- see tryFlipCoplanarQuadDiagonal()'s
    // own comment for exactly which cases phase 1 does and does not cover,
    // and this class's header for why full generality past a genuinely
    // non-convex (Schonhardt-style) boundary remains open research.
    FacetRecoveryResult recoverFacets(const std::vector<std::array<std::size_t, 3>>& facets,
                                       int maxRounds = 4);

    // Carves an interior cavity out of an already-tetrahedralized domain:
    // removes every tetrahedron reachable from the tetrahedron containing
    // `seed`, walking face-adjacency, without crossing any face present in
    // `walls` (typically a FacetRecoveryResult's own recoveredFacets) --
    // the standard "hole point" technique for turning a solid
    // tetrahedralization of an outer-plus-inner boundary into one that
    // actually respects the inner boundary as a void. Returns the number
    // of tetrahedra removed (0 if `seed` isn't inside any tetrahedron).
    std::size_t removeRegion(const core::Vector3& seed, const std::vector<std::array<std::size_t, 3>>& walls);

    std::size_t pointCount() const { return points_.size(); }
    std::size_t tetrahedronCount() const { return tetrahedra_.size(); }
    const core::Vector3& point(std::size_t i) const { return points_.at(i); }
    const Tetrahedron& tetrahedron(std::size_t i) const { return tetrahedra_.at(i); }

    // Brute-force check of the defining Delaunay property: no point lies
    // strictly inside any tetrahedron's circumsphere. O(n * tetrahedronCount);
    // intended for validating modest point sets in tests, not as a
    // runtime query.
    bool satisfiesDelaunayProperty(double tolerance = 1e-6) const;

private:
    double inCircumsphere(const Tetrahedron& t, std::size_t pointIndex) const;

    // True iff `face` (any vertex order) is a face of some tetrahedron.
    bool hasFace(const std::array<std::size_t, 3>& face) const;

    // recoverFacets()'s per-facet recursion: checks `facet`, and if
    // missing and roundsLeft > 0, inserts its centroid and recurses on the
    // three resulting sub-triangles.
    void recoverFacetRecursive(const std::array<std::size_t, 3>& facet, int roundsLeft,
                                std::vector<std::array<std::size_t, 3>>& recovered,
                                std::vector<std::array<std::size_t, 3>>& unrecovered);

    // recoverFacets()'s combinatorial first phase for one candidate pair.
    // Applicable when `f1` and `f2` share exactly one edge and, together
    // with that edge's two endpoints, describe a single coplanar
    // quadrilateral -- i.e. f1/f2 are the two triangles of one diagonal
    // choice on a quad whose *other* diagonal is what the tetrahedralization
    // actually built. Recovers both at once by replacing the two
    // tetrahedra straddling that quad (which must share one common fourth
    // vertex -- the simple "bipyramid" case) with two others built on the
    // requested diagonal instead, after confirming explicitly (not by
    // trusting makeTetrahedron()'s own orientation auto-correction) that
    // the quad is convex as seen from that shared vertex, so the swap
    // covers exactly the same volume rather than overlapping or leaving a
    // gap.
    //
    // Deliberately narrow: a facet pair not fitting this exact shape (no
    // shared edge, not coplanar, the two straddling tetrahedra having
    // *different* fourth vertices, or a non-convex quad) is left for
    // recoverFacetRecursive()'s Steiner fallback rather than guessed at --
    // the general case (an arbitrary local cavity needing more than two
    // tetrahedra retetrahedralized) is the genuinely research-level part of
    // this problem this class's header already flags.
    //
    // Returns true and mutates the mesh iff the flip was performed.
    bool tryFlipCoplanarQuadDiagonal(const std::array<std::size_t, 3>& f1, const std::array<std::size_t, 3>& f2);

    // Shared cavity-finding step between tetrahedralize()'s main loop and
    // insertSteinerPoint(): every tetrahedron whose circumsphere contains
    // point `p`, and the triangular faces bounding the star-shaped cavity
    // they leave behind (each face belongs to exactly one bad tetrahedron).
    //
    // considerHullVisibility enables the extra hull-face/ghost-tetrahedron
    // logic insertSteinerPoint() needs to support points outside the
    // current hull (see its own comment). tetrahedralize()'s own loop
    // always passes false: every point it inserts is strictly interior to
    // the enclosing super-tetrahedron by construction, so that logic is
    // never geometrically needed there, and false is the simpler,
    // longer-tested code path.
    void findCavity(std::size_t p, bool considerHullVisibility, std::vector<std::size_t>& badIndices,
                     std::vector<std::array<std::size_t, 3>>& boundary) const;

    std::vector<core::Vector3> points_;
    std::vector<Tetrahedron> tetrahedra_;
};

} // namespace aether::mesh
