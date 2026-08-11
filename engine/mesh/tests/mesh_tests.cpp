#include "aether/mesh/DelaunayTetrahedralization3D.hpp"
#include "aether/mesh/DelaunayTriangulation2D.hpp"
#include "aether/mesh/PolygonTriangulation2D.hpp"
#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>

using namespace aether::core;
using namespace aether::mesh;

namespace {

bool nearlyEqual(double a, double b, double tolerance = 1e-9) {
    return std::fabs(a - b) <= tolerance;
}

void testBasicLayout() {
    StructuredGrid3D grid(Vector3(0, 0, 0), Vector3(2, 1, 1), 4, 2, 1);
    AETHER_CHECK(grid.cellCount() == 8);
    AETHER_CHECK(nearlyEqual(grid.spacing().x, 0.5));
    AETHER_CHECK(nearlyEqual(grid.spacing().y, 0.5));
    AETHER_CHECK(nearlyEqual(grid.spacing().z, 1.0));
    AETHER_CHECK(nearlyEqual(grid.cellVolume(), 0.5 * 0.5 * 1.0));

    const Vector3 center000 = grid.cellCenter(0, 0, 0);
    AETHER_CHECK(nearlyEqual(center000.x, 0.25));
    AETHER_CHECK(nearlyEqual(center000.y, 0.25));
    AETHER_CHECK(nearlyEqual(center000.z, 0.5));
}

void testCellIndexIsUnique() {
    StructuredGrid3D grid(Vector3(0, 0, 0), Vector3(1, 1, 1), 3, 3, 3);
    bool seen[27] = {false};
    for (std::size_t k = 0; k < 3; ++k) {
        for (std::size_t j = 0; j < 3; ++j) {
            for (std::size_t i = 0; i < 3; ++i) {
                const std::size_t idx = grid.cellIndex(i, j, k);
                AETHER_CHECK(idx < 27);
                AETHER_CHECK(!seen[idx]);
                seen[idx] = true;
            }
        }
    }
}

void testNeighborBounds() {
    StructuredGrid3D grid(Vector3(0, 0, 0), Vector3(1, 1, 1), 3, 3, 3);
    AETHER_CHECK(!grid.hasNeighbor(0, 0, 0, -1, 0, 0));
    AETHER_CHECK(grid.hasNeighbor(0, 0, 0, 1, 0, 0));
    AETHER_CHECK(!grid.hasNeighbor(2, 2, 2, 1, 0, 0));
    AETHER_CHECK(grid.hasNeighbor(1, 1, 1, -1, 0, 0));
}

// A small point set in "general position" (irregular, deterministic
// jitter added to a grid so no four points are exactly co-circular --
// a perfectly regular grid is a classic degenerate case for Delaunay
// triangulation, with ties/numerical fragility that would make this a
// worse first test, not a better one).
void testDelaunayTriangulationSatisfiesDelaunayProperty() {
    DelaunayTriangulation2D triangulation;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < 5; ++j) {
            const double x = static_cast<double>(i) + 0.1 * std::sin(3.7 * i + 1.3 * j);
            const double y = static_cast<double>(j) + 0.1 * std::cos(2.1 * i + 4.2 * j);
            triangulation.addPoint(x, y);
        }
    }
    triangulation.triangulate();

    AETHER_CHECK(triangulation.pointCount() == 25);
    AETHER_CHECK(triangulation.triangleCount() > 0);
    AETHER_CHECK(triangulation.satisfiesDelaunayProperty());

    // No degenerate (zero-area) triangles, and every triangle CCW-wound.
    for (std::size_t t = 0; t < triangulation.triangleCount(); ++t) {
        const auto& tri = triangulation.triangle(t);
        const Vector3& a = triangulation.point(tri.vertices[0]);
        const Vector3& b = triangulation.point(tri.vertices[1]);
        const Vector3& c = triangulation.point(tri.vertices[2]);
        const double signedArea2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        AETHER_CHECK(signedArea2 > 1e-9);
    }
}

// Regression test for the too-small-super-triangle bug -- the 2D twin of
// testTetrahedralizationFillsConvexHullOnHullAdjacentSliverCase below.
// The super triangle's vertices are a *finite* stand-in for points at
// infinity; when they sit too close, a triangle built on real points plus
// a super vertex gets a finite circumcircle that can wrongly test as
// containing a nearby real point, deleting a triangle that should have
// survived. See triangulate()'s own comment for the full mechanism and
// the A/B measurement (the old 20*deltaMax scale failed 10 of 200
// randomized point sets; the current scale passes all 200).
//
// These 14 points are one of those measured failures, kept verbatim. The
// invariant checked is coverage, not just validity: a triangulation of a
// point set must exactly fill its convex hull, so the summed triangle
// area must equal the hull area (computed independently via qhull's
// ConvexHull, which uses no Delaunay algorithm at all). With the old
// scale this produced 18 triangles totalling 11.184472..., short of the
// hull by one triangle's worth of area.
//
// Note satisfiesDelaunayProperty() does NOT catch this: the buggy output
// satisfied it perfectly. Circumcircle-emptiness says nothing about
// whether the triangles actually cover the hull.
void testTriangulationFillsConvexHullOnHullAdjacentSliverCase() {
    DelaunayTriangulation2D triangulation;
    triangulation.addPoint(-1.3188690253489734, 1.4770029849644346);
    triangulation.addPoint(-1.689162684159931, 1.903003372933219);
    triangulation.addPoint(-0.9285840630165665, 1.2221582906384008);
    triangulation.addPoint(0.011479057953599714, 0.7869403196580191);
    triangulation.addPoint(1.766025650841355, 0.23036728694670439);
    triangulation.addPoint(1.0400022097229087, -1.7478773565584969);
    triangulation.addPoint(-1.900102265528079, 1.423962424326659);
    triangulation.addPoint(1.3633210999009395, 1.8851982382895422);
    triangulation.addPoint(1.4816110656953922, -1.1435436331524156);
    triangulation.addPoint(1.2143882434656375, 0.5702636158209788);
    triangulation.addPoint(0.35214282083424164, -0.9174072883030893);
    triangulation.addPoint(-0.9836596971421101, -1.8674129393947223);
    triangulation.addPoint(0.706064541872931, -0.2728510605701189);
    triangulation.addPoint(1.4847135912224152, 1.3681442520725944);
    triangulation.triangulate();

    AETHER_CHECK(triangulation.satisfiesDelaunayProperty());
    AETHER_CHECK(triangulation.triangleCount() == 19); // was 18 before the fix

    double totalArea = 0.0;
    for (std::size_t t = 0; t < triangulation.triangleCount(); ++t) {
        const auto& tri = triangulation.triangle(t);
        const Vector3& a = triangulation.point(tri.vertices[0]);
        const Vector3& b = triangulation.point(tri.vertices[1]);
        const Vector3& c = triangulation.point(tri.vertices[2]);
        totalArea += ((b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)) / 2.0;
    }
    AETHER_CHECK(nearlyEqual(totalArea, 11.188140298284395, 1e-12));
}

// A convex point set (regular pentagon plus its center) has a known
// triangle count: triangulating n points with h of them on the convex
// hull always yields 2n - 2 - h triangles, regardless of which valid
// Delaunay configuration is chosen among any co-circular ties. Here
// n=6, h=5 (the pentagon; the center point is interior), so exactly
// 2*6-2-5=5 triangles are expected -- a second, independent structural
// check beyond the circumcircle property itself.
void testDelaunayTriangulationOfConvexPolygonPlusCenterMatchesEulerCount() {
    const double kPi = 3.14159265358979323846;
    DelaunayTriangulation2D triangulation;
    for (int i = 0; i < 5; ++i) {
        const double angle = 2.0 * kPi * static_cast<double>(i) / 5.0;
        triangulation.addPoint(std::cos(angle), std::sin(angle));
    }
    triangulation.addPoint(0.0, 0.0);
    triangulation.triangulate();

    AETHER_CHECK(triangulation.satisfiesDelaunayProperty());
    AETHER_CHECK(triangulation.triangleCount() == 5);
}

// L-shaped hexagon (concave at one vertex) -- a convex polygon wouldn't
// exercise ear clipping's reflex-vertex handling at all, so this is the
// meaningful first test, not just a convenient one. Validated two exact,
// self-contained ways: (1) triangle count for any simple n-vertex polygon
// is always n-2, here 6-2=4; (2) the sum of triangle areas must equal the
// polygon's own shoelace-formula area exactly -- if ear clipping ever
// produced overlapping or gapped triangles, this sum would be wrong.
void testPolygonTriangulationOfConcavePolygonPartitionsExactly() {
    PolygonTriangulation2D triangulation;
    triangulation.addVertex(0.0, 0.0);
    triangulation.addVertex(2.0, 0.0);
    triangulation.addVertex(2.0, 1.0);
    triangulation.addVertex(1.0, 1.0);
    triangulation.addVertex(1.0, 2.0);
    triangulation.addVertex(0.0, 2.0);
    triangulation.triangulate();

    AETHER_CHECK(triangulation.triangleCount() == 4); // n - 2 for a simple polygon
    AETHER_CHECK(triangulation.isLocallyDelaunay());

    const double expectedArea = triangulation.polygonArea();
    AETHER_CHECK(nearlyEqual(expectedArea, 3.0, 1e-9)); // 2x1 rectangle + 1x1 square

    double summedTriangleArea = 0.0;
    for (std::size_t t = 0; t < triangulation.triangleCount(); ++t) {
        const auto& tri = triangulation.triangle(t);
        const Vector3& a = triangulation.vertex(tri.vertices[0]);
        const Vector3& b = triangulation.vertex(tri.vertices[1]);
        const Vector3& c = triangulation.vertex(tri.vertices[2]);
        const double signedArea2 = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
        AETHER_CHECK(signedArea2 > 1e-9); // non-degenerate and CCW
        summedTriangleArea += signedArea2 / 2.0;
    }
    AETHER_CHECK(nearlyEqual(summedTriangleArea, expectedArea, 1e-9));
}

// A convex quadrilateral A=(0,0), B=(3,0), C=(3,1), D=(0,4) chosen because
// plain ear clipping (tried by increasing vertex index) picks vertex A as
// the first ear, which forces the diagonal B-D -- but that diagonal is
// *not* locally Delaunay here (confirmed empirically before writing this
// test), so triangulate()'s flip pass must replace it with A-C. This is a
// concrete demonstration that the flip step actually changes the result,
// not just a property that happens to already hold.
void testPolygonTriangulationFlipsToLocallyDelaunayDiagonal() {
    PolygonTriangulation2D triangulation;
    triangulation.addVertex(0.0, 0.0); // A = 0
    triangulation.addVertex(3.0, 0.0); // B = 1
    triangulation.addVertex(3.0, 1.0); // C = 2
    triangulation.addVertex(0.0, 4.0); // D = 3
    triangulation.triangulate();

    AETHER_CHECK(triangulation.triangleCount() == 2);
    AETHER_CHECK(triangulation.isLocallyDelaunay());

    bool usesDiagonalAC = false;
    for (std::size_t t = 0; t < triangulation.triangleCount(); ++t) {
        const auto& v = triangulation.triangle(t).vertices;
        const bool hasA = v[0] == 0 || v[1] == 0 || v[2] == 0;
        const bool hasC = v[0] == 2 || v[1] == 2 || v[2] == 2;
        if (hasA && hasC) {
            usesDiagonalAC = true;
        }
    }
    AETHER_CHECK(usesDiagonalAC); // confirms the flip actually happened, not a no-op
}

// A small 3D point set in "general position" (jittered 3x3x3 grid), the
// direct 3D analog of the 5x5 jittered-grid 2D test above -- a perfectly
// regular grid would put many points exactly co-spherical, a degenerate
// case that's a worse first test than a better one.
void testDelaunayTetrahedralizationSatisfiesDelaunayProperty() {
    DelaunayTetrahedralization3D tet;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                const double x = static_cast<double>(i) + 0.1 * std::sin(3.7 * i + 1.3 * j + 2.4 * k);
                const double y = static_cast<double>(j) + 0.1 * std::cos(2.1 * i + 4.2 * j + 1.1 * k);
                const double z = static_cast<double>(k) + 0.1 * std::sin(1.6 * i + 2.9 * j + 3.3 * k);
                tet.addPoint(x, y, z);
            }
        }
    }
    tet.tetrahedralize();

    AETHER_CHECK(tet.pointCount() == 27);
    AETHER_CHECK(tet.tetrahedronCount() > 0);
    AETHER_CHECK(tet.satisfiesDelaunayProperty());

    // No degenerate (zero-volume) tetrahedra, and every one positively
    // oriented (signed volume > 0, as tetrahedralize() guarantees).
    for (std::size_t t = 0; t < tet.tetrahedronCount(); ++t) {
        const auto& tv = tet.tetrahedron(t);
        const Vector3& a = tet.point(tv.vertices[0]);
        const Vector3& b = tet.point(tv.vertices[1]);
        const Vector3& c = tet.point(tv.vertices[2]);
        const Vector3& d = tet.point(tv.vertices[3]);
        const double signedVolume6 = (b - a).dot((c - a).cross(d - a));
        AETHER_CHECK(signedVolume6 > 1e-9);
    }
}

// A convex hexahedron with an exactly known volume: a square bipyramid
// (base vertices at (+-1,0,0) and (0,+-1,0), apexes at (0,0,2) and
// (0,0,-1)). Unequal apex heights (and both different from the base's
// own radius) avoid the accidental extra co-sphericity that a regular
// octahedron would have. Volume = (1/3)*baseArea*(h_top + h_bottom):
// the base is a square with diagonals of length 2, area = 2*2/2 = 2, so
// volume = (1/3)*2*(2+1) = 2. If tetrahedralize() ever produced
// overlapping or gapped tetrahedra, the summed |signed volume| would not
// match this exactly -- the 3D analog of the 2D "summed triangle area
// equals shoelace area" invariant.
void testDelaunayTetrahedralizationOfBipyramidPartitionsExactVolume() {
    DelaunayTetrahedralization3D tet;
    tet.addPoint(1.0, 0.0, 0.0);
    tet.addPoint(0.0, 1.0, 0.0);
    tet.addPoint(-1.0, 0.0, 0.0);
    tet.addPoint(0.0, -1.0, 0.0);
    tet.addPoint(0.0, 0.0, 2.0);
    tet.addPoint(0.0, 0.0, -1.0);
    tet.tetrahedralize();

    AETHER_CHECK(tet.satisfiesDelaunayProperty());
    AETHER_CHECK(tet.tetrahedronCount() > 0);

    double summedVolume = 0.0;
    for (std::size_t t = 0; t < tet.tetrahedronCount(); ++t) {
        const auto& tv = tet.tetrahedron(t);
        const Vector3& a = tet.point(tv.vertices[0]);
        const Vector3& b = tet.point(tv.vertices[1]);
        const Vector3& c = tet.point(tv.vertices[2]);
        const Vector3& d = tet.point(tv.vertices[3]);
        const double signedVolume6 = (b - a).dot((c - a).cross(d - a));
        AETHER_CHECK(signedVolume6 > 1e-9); // non-degenerate and positively oriented
        summedVolume += signedVolume6 / 6.0;
    }
    AETHER_CHECK(nearlyEqual(summedVolume, 2.0, 1e-9));
}

double sumTetrahedraVolume(const DelaunayTetrahedralization3D& tet) {
    double total = 0.0;
    for (std::size_t t = 0; t < tet.tetrahedronCount(); ++t) {
        const auto& tv = tet.tetrahedron(t);
        const Vector3& a = tet.point(tv.vertices[0]);
        const Vector3& b = tet.point(tv.vertices[1]);
        const Vector3& c = tet.point(tv.vertices[2]);
        const Vector3& d = tet.point(tv.vertices[3]);
        total += (b - a).dot((c - a).cross(d - a)) / 6.0;
    }
    return total;
}

// insertSteinerPoint() for an *interior* point validated the same way
// tetrahedralize() itself is: the Delaunay property must still hold
// afterward, and -- since inserting an interior point only subdivides
// existing volume, never adds or removes any -- the summed tetrahedra
// volume must stay exactly the same (the same "measure an exact
// invariant" style as the bipyramid test above, now applied across an
// incremental update instead of a single full build).
void testSteinerPointInsertionConservesVolumeAndDelaunayProperty() {
    DelaunayTetrahedralization3D tet;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                const double x = static_cast<double>(i) + 0.1 * std::sin(3.7 * i + 1.3 * j + 2.4 * k);
                const double y = static_cast<double>(j) + 0.1 * std::cos(2.1 * i + 4.2 * j + 1.1 * k);
                const double z = static_cast<double>(k) + 0.1 * std::sin(1.6 * i + 2.9 * j + 3.3 * k);
                tet.addPoint(x, y, z);
            }
        }
    }
    tet.tetrahedralize();
    AETHER_CHECK(tet.satisfiesDelaunayProperty());

    const double volumeBefore = sumTetrahedraVolume(tet);
    const std::size_t tetCountBefore = tet.tetrahedronCount();
    const std::size_t pointCountBefore = tet.pointCount();

    // The grid's own centroid: well inside the jittered 3x3x3 cloud's hull.
    AETHER_CHECK(tet.insertSteinerPoint(1.0, 1.0, 1.0));
    AETHER_CHECK(tet.pointCount() == pointCountBefore + 1);
    AETHER_CHECK(tet.tetrahedronCount() > tetCountBefore); // net subdivision, strictly more tets
    AETHER_CHECK(tet.satisfiesDelaunayProperty());
    AETHER_CHECK(nearlyEqual(sumTetrahedraVolume(tet), volumeBefore, 1e-9));
}

// insertSteinerPoint() for a point *outside* the current hull: the exact
// opposite volume invariant from the interior case above -- extending the
// hull necessarily *adds* volume, it doesn't just subdivide what's already
// there -- checked here against an exactly known number rather than just
// "some increase happened". The right-angle tetrahedron O=(0,0,0),
// X=(1,0,0), Y=(0,1,0), Z=(0,0,1) has volume exactly 1/6; inserting
// p=(1,1,1) (straight out through the face XYZ, the one opposite the
// origin) must leave that original tetrahedron completely untouched (its
// own circumsphere doesn't contain p -- confirmed independently in Python
// during development, center ~(0.5,0.5,0.5)-ish, p well outside it) and
// add exactly one new cap tetrahedron (X,Y,Z,p), whose volume is exactly
// 1/3 (computed the same way: half the volume of the unit cube corner
// pyramid XYZp). Total after insertion: exactly 1/6 + 1/3 = 1/2.
//
// This exact case is also what caught a real bug during development: an
// earlier version of findCavity() marked the *entire* original tetrahedron
// bad whenever it owned a face visible from p, discarding it and re-fanning
// all of its non-visible faces to p -- a plausible-looking but wrong rule
// (it still produced watertight, correctly-summed-volume output, so the
// bug survived a pure volume check!) that was only caught by cross-checking
// against scipy's own Delaunay/ConvexHull implementation on this exact
// point set, which showed the original tetrahedron must stay untouched
// with only one new tetrahedron added alongside it -- exactly what this
// test asserts, plus the Delaunay property, which is what the old version
// actually violated.
void testSteinerPointInsertionExtendsHullWithExactVolume() {
    DelaunayTetrahedralization3D tet;
    tet.addPoint(0.0, 0.0, 0.0);
    tet.addPoint(1.0, 0.0, 0.0);
    tet.addPoint(0.0, 1.0, 0.0);
    tet.addPoint(0.0, 0.0, 1.0);
    tet.tetrahedralize();
    AETHER_CHECK(tet.satisfiesDelaunayProperty());
    AETHER_CHECK(nearlyEqual(sumTetrahedraVolume(tet), 1.0 / 6.0, 1e-12));

    AETHER_CHECK(tet.insertSteinerPoint(1.0, 1.0, 1.0));
    AETHER_CHECK(tet.pointCount() == 5);
    AETHER_CHECK(tet.tetrahedronCount() == 2); // original tet kept, one new cap added
    AETHER_CHECK(tet.satisfiesDelaunayProperty());
    AETHER_CHECK(nearlyEqual(sumTetrahedraVolume(tet), 0.5, 1e-12));
}

// Regression test for the "missing hull-adjacent tetrahedron" bug fixed by
// enlarging the super-tetrahedron (see tetrahedralize()'s own comment for
// the mechanism, and the class header for how the first diagnosis of this
// was wrong). These are the exact 7 points that first exposed it, kept
// verbatim: with the old too-small super-tetrahedron this produced 6
// tetrahedra totalling 2.708055... , short of the true convex hull volume
// by the volume of one thin, valid, hull-adjacent tetrahedron.
//
// The expected total is the convex hull volume, since a tetrahedralization
// of a point set must exactly fill its convex hull -- checked here against
// the independently-computed value (2.709012828039542, from qhull via
// scipy.spatial.ConvexHull, which needs no Delaunay algorithm at all and
// so is a genuinely independent reference rather than a self-comparison).
// Note satisfiesDelaunayProperty() alone would NOT catch this: the buggy
// output satisfied it perfectly, because circumsphere-emptiness says
// nothing about whether the tetrahedra actually cover the hull.
void testTetrahedralizationFillsConvexHullOnHullAdjacentSliverCase() {
    DelaunayTetrahedralization3D tet;
    tet.addPoint(1.4737617982207118, -0.8528567809098169, -0.09741757121024008);
    tet.addPoint(0.8818641874977176, -1.5726403556472421, -1.5119823378293402);
    tet.addPoint(1.7058871049329896, -0.5884653761446414, 0.4584386055893468);
    tet.addPoint(-0.21987487313532084, -0.5879755232108659, -0.548828267652588);
    tet.addPoint(-0.4376550633723242, -1.8110581604919749, 0.034080232074870764);
    tet.addPoint(1.540759526282839, -1.6205730261557352, 0.6577405728359289);
    tet.addPoint(0.9445458377233287, 1.12204139437512, -1.1761955929033237);
    tet.tetrahedralize();

    AETHER_CHECK(tet.satisfiesDelaunayProperty());
    AETHER_CHECK(tet.tetrahedronCount() == 7); // was 6 before the fix
    AETHER_CHECK(nearlyEqual(sumTetrahedraVolume(tet), 2.709012828039542, 1e-12));
}

// Builds the 6 vertices of a regular octahedron ("cross-polytope")
// centered at `center` with axis-vertex distance `r`: index 0/1 = +x/-x,
// 2/3 = +y/-y, 4/5 = +z/-z (relative to center).
void addOctahedronPoints(DelaunayTetrahedralization3D& tet, const Vector3& center, double r) {
    tet.addPoint(center.x + r, center.y, center.z);
    tet.addPoint(center.x - r, center.y, center.z);
    tet.addPoint(center.x, center.y + r, center.z);
    tet.addPoint(center.x, center.y - r, center.z);
    tet.addPoint(center.x, center.y, center.z + r);
    tet.addPoint(center.x, center.y, center.z - r);
}

// The 8 triangular faces of an octahedron built by addOctahedronPoints(),
// vertex indices offset by `base`: every combination of one vertex per
// axis. Deliberately triangles, not quads split by an arbitrary diagonal
// choice -- see the class header and
// testTetrahedralization3DHonestlyReportsUnrecoverableCoplanarQuadFacets()
// for why that distinction is the whole point of this test's shape choice.
std::vector<std::array<std::size_t, 3>> octahedronFaces(std::size_t base) {
    std::vector<std::array<std::size_t, 3>> faces;
    for (std::size_t x : {base + 0, base + 1}) {
        for (std::size_t y : {base + 2, base + 3}) {
            for (std::size_t z : {base + 4, base + 5}) {
                faces.push_back({x, y, z});
            }
        }
    }
    return faces;
}

// **The end-to-end validation of task #74's tractable slice**: an outer
// octahedron (R=2) with a smaller octahedron-shaped hole (r=0.5) carved out
// of its middle, both centered at the same point. Octahedron faces are
// genuine triangles (no quad-diagonal choice), so this shape sidesteps the
// coplanar-facet ambiguity the sibling test below documents, and lets the
// full recover-then-carve pipeline be checked against an exact ground
// truth: the L1-ball volume formula (4/3)*R^3 for an octahedron with
// axis-vertex distance R, so the expected carved volume is
// (4/3)*(R^3 - r^3) = (4/3)*(8 - 0.125) = 10.5 exactly.
//
// Measured directly before writing this test (not guessed): with this
// shape, missingFacets() already returns 0 of the 16 requested facets --
// the plain unconstrained tetrahedralization recovers every octahedron
// face on its own, since a triangular face never needs the kind of local
// re-triangulation a planar quad does. recoverFacets() is still called (it
// is a no-op here, confirmed by recoveredFacets covering exactly those 16
// and unrecovered being empty) so the full pipeline -- not just the lucky
// case -- is what is under test.
void testTetrahedralization3DRecoversHollowOctahedronFacetsAndCarvesHoleWithExactVolume() {
    DelaunayTetrahedralization3D tet;
    const Vector3 center(1.0, 1.0, 1.0);
    addOctahedronPoints(tet, center, 2.0);
    addOctahedronPoints(tet, center, 0.5);
    tet.tetrahedralize();
    AETHER_CHECK(tet.satisfiesDelaunayProperty());

    std::vector<std::array<std::size_t, 3>> facets = octahedronFaces(0);
    const std::vector<std::array<std::size_t, 3>> innerFacets = octahedronFaces(6);
    facets.insert(facets.end(), innerFacets.begin(), innerFacets.end());
    AETHER_CHECK(facets.size() == 16);

    AETHER_CHECK(tet.missingFacets(facets).empty()); // measured: recovered automatically

    const auto result = tet.recoverFacets(facets, 4);
    AETHER_CHECK(result.recoveredFacets.size() == 16);
    AETHER_CHECK(result.unrecovered.empty());

    const std::size_t removed = tet.removeRegion(center, result.recoveredFacets);
    AETHER_CHECK(removed > 0);
    AETHER_CHECK(tet.satisfiesDelaunayProperty()); // removal only deletes tetrahedra, never adds points

    const double expectedVolume = (4.0 / 3.0) * (2.0 * 2.0 * 2.0 - 0.5 * 0.5 * 0.5);
    AETHER_CHECK(nearlyEqual(sumTetrahedraVolume(tet), expectedVolume, 1e-9)); // measured: exactly 10.5
}

// **An honest, documented limitation, found by testing, not hidden.** A
// single axis-aligned cube's 6 square faces, each split into 2 triangles
// by picking one diagonal: the plain tetrahedralization of just the 8
// corners is free to pick *either* diagonal for each face (both give an
// equally valid, empty-circumsphere local triangulation of that exactly
// planar 4-point set -- the same kind of tie DelaunayTriangulation2D's own
// jittered-grid test deliberately avoids by jittering, except here the
// points genuinely must stay exactly coplanar to be a real cube face).
// Measured directly: this specific cube ends up missing exactly 6 of the
// 12 requested half-square triangles (the other diagonal was chosen for
// those 3 faces).
//
// recoverFacets()'s centroid-Steiner-point heuristic does not recover
// these: a missing facet's centroid lies exactly on the same plane as its
// face's 4 corners (by construction, since the facet is coplanar with its
// diagonal partner), so inserting it either gets rejected immediately by
// insertSteinerPoint()'s own degenerate-orientation guard, or (as observed
// during development on a related hollow-cube case) recurses through
// several further exactly-coplanar sub-triangles before still failing --
// this project's practice is to report that honestly via `unrecovered`
// rather than silently leaving the caller with an incomplete boundary, so
// that is exactly what is checked here: nothing is lost or corrupted, both
// lists just partition the 12 input facets deterministically.
void testTetrahedralization3DHonestlyReportsUnrecoverableCoplanarQuadFacets() {
    DelaunayTetrahedralization3D tet;
    tet.addPoint(0.0, 0.0, 0.0);
    tet.addPoint(2.0, 0.0, 0.0);
    tet.addPoint(2.0, 2.0, 0.0);
    tet.addPoint(0.0, 2.0, 0.0);
    tet.addPoint(0.0, 0.0, 2.0);
    tet.addPoint(2.0, 0.0, 2.0);
    tet.addPoint(2.0, 2.0, 2.0);
    tet.addPoint(0.0, 2.0, 2.0);
    tet.tetrahedralize();

    const std::vector<std::array<std::size_t, 3>> facets = {
        {0, 1, 2}, {0, 2, 3}, {4, 5, 6}, {4, 6, 7}, {0, 1, 5}, {0, 5, 4},
        {3, 2, 6}, {3, 6, 7}, {0, 3, 7}, {0, 7, 4}, {1, 2, 6}, {1, 6, 5},
    };

    AETHER_CHECK(tet.missingFacets(facets).size() == 6); // measured directly

    const auto result = tet.recoverFacets(facets, 4);
    // Every input facet is accounted for exactly once, in one list or the
    // other -- nothing silently dropped.
    AETHER_CHECK(result.recoveredFacets.size() + result.unrecovered.size() == facets.size());
    AETHER_CHECK(result.recoveredFacets.size() == 6); // the 6 already present
    AETHER_CHECK(result.unrecovered.size() == 6);      // measured: the heuristic never recovers these
}

} // namespace

int main() {
    testBasicLayout();
    testCellIndexIsUnique();
    testNeighborBounds();
    testDelaunayTriangulationSatisfiesDelaunayProperty();
    testTriangulationFillsConvexHullOnHullAdjacentSliverCase();
    testDelaunayTriangulationOfConvexPolygonPlusCenterMatchesEulerCount();
    testPolygonTriangulationOfConcavePolygonPartitionsExactly();
    testPolygonTriangulationFlipsToLocallyDelaunayDiagonal();
    testDelaunayTetrahedralizationSatisfiesDelaunayProperty();
    testDelaunayTetrahedralizationOfBipyramidPartitionsExactVolume();
    testSteinerPointInsertionConservesVolumeAndDelaunayProperty();
    testSteinerPointInsertionExtendsHullWithExactVolume();
    testTetrahedralizationFillsConvexHullOnHullAdjacentSliverCase();
    testTetrahedralization3DRecoversHollowOctahedronFacetsAndCarvesHoleWithExactVolume();
    testTetrahedralization3DHonestlyReportsUnrecoverableCoplanarQuadFacets();
    std::puts("aether_mesh_tests: all tests passed");
    return 0;
}
