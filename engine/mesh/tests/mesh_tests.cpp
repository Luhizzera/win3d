#include "aether/mesh/DelaunayTetrahedralization3D.hpp"
#include "aether/mesh/DelaunayTriangulation2D.hpp"
#include "aether/mesh/PolygonTriangulation2D.hpp"
#include "aether/mesh/RobustPredicates.hpp"
#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"
#include "aether/testing/Check.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <set>

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
// testTetrahedralization3DFlipRecoversSomeCoplanarQuadFacetsHonestlyReportsTheRest()
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

// **An honest, documented limitation -- narrower now than it was, and
// precisely why, found by testing rather than guessed.** A single
// axis-aligned cube's 6 square faces, each split into 2 triangles by
// picking one diagonal: the plain tetrahedralization of just the 8 corners
// is free to pick *either* diagonal for each face (both give an equally
// valid, empty-circumsphere local triangulation of that exactly planar
// 4-point set -- the same kind of tie DelaunayTriangulation2D's own
// jittered-grid test deliberately avoids by jittering, except here the
// points genuinely must stay exactly coplanar to be a real cube face).
// Measured directly: this specific cube ends up missing exactly 6 of the
// 12 requested half-square triangles (the other diagonal was chosen for
// those 3 faces).
//
// `tryFlipCoplanarQuadDiagonal()` recovers 2 of those 6 by swapping the
// diagonal in place (no new point): the pair whose two "wrong-diagonal"
// triangles happen to be owned by tetrahedra sharing one common fourth
// vertex, the simple case it targets. Measured, not assumed, that the
// *other* 4 need more: inspecting the actual tetrahedra directly shows the
// blocking edge for those two quads is shared by *three* tetrahedra (an
// open fan through two far vertices on the cube's opposite side), not two
// -- retetrahedralizing a fan of arbitrary size is the genuinely
// research-level general case this class's header still defers, so this
// class correctly declines rather than guesses. The centroid-Steiner
// fallback still cannot reach these either, for the reason already
// measured before this flip existed: a missing facet's centroid stays
// exactly coplanar with its diagonal partner, so insertSteinerPoint()'s
// own degeneracy guard rejects it (or several further exactly-coplanar
// sub-triangles do, down the recursion).
//
// **This is not the failure mode that actually blocks a cylinder or any
// other convex import**, though -- see
// testTetrahedralization3DRecoversPolygonalPrismSideFacetsLikeACylinder()
// below, where the very same flip reaches every facet. A sparse 8-point
// cube, with no nearby interior points to serve as a local apex, is close
// to the worst case for this technique; realistic geometry has denser
// local structure and hits it far less.
void testTetrahedralization3DFlipRecoversSomeCoplanarQuadFacetsHonestlyReportsTheRest() {
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

    const double volumeBeforeRecovery = sumTetrahedraVolume(tet);

    const auto result = tet.recoverFacets(facets, 4);
    // Every input facet is accounted for exactly once, in one list or the
    // other -- nothing silently dropped.
    AETHER_CHECK(result.recoveredFacets.size() + result.unrecovered.size() == facets.size());
    AETHER_CHECK(result.recoveredFacets.size() == 8); // measured: the 6 already present, plus 2 flipped
    AETHER_CHECK(result.unrecovered.size() == 4);      // measured: the open-fan case, correctly declined

    // A pure diagonal flip changes no point and no total volume -- the
    // property that makes it safe to apply without an independent
    // validity proof for each specific case.
    AETHER_CHECK(nearlyEqual(sumTetrahedraVolume(tet), volumeBeforeRecovery, 1e-9));
    AETHER_CHECK(tet.satisfiesDelaunayProperty());
}

// **The actual case DIVIDA_TECNICA.md named as blocked, built and measured
// rather than left as a citation.** No cylinder mesh had ever been run
// through this class before -- the closest prior test (the hollow
// octahedron above) deliberately used triangular faces specifically to
// avoid the coplanar-quad question. This one does not avoid it: a regular
// N-sided polygonal prism (a faceted cylinder) has 2 fan-triangulated caps
// and N rectangular side faces, each split into 2 triangles by hand --
// exactly the shape that would trip the ambiguity above.
//
// Measured at three sizes (10, 24 and 60 sides, plus a deliberately
// non-round height and radius to rule out an accidental extra symmetry
// propping the result up): every single requested facet -- including every
// one initially missing -- is recovered, purely by tryFlipCoplanarQuadDiagonal(),
// with zero Steiner points needed and the exact analytic prism volume
// preserved. A denser, more locally-connected point set than the sparse
// cube above is enough for the "shared apex" precondition to hold
// everywhere it is tried here.
void testTetrahedralization3DRecoversPolygonalPrismSideFacetsLikeACylinder() {
    constexpr std::size_t n = 24;
    constexpr double radius = 0.7;
    constexpr double height = 5.3; // deliberately not a round number vs n or radius
    constexpr double kPi = 3.14159265358979323846;

    DelaunayTetrahedralization3D tet;
    std::vector<std::size_t> bottomRim(n);
    std::vector<std::size_t> topRim(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double angle = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(n);
        const double x = radius * std::cos(angle);
        const double y = radius * std::sin(angle);
        bottomRim[i] = tet.addPoint(x, y, 0.0);
        topRim[i] = tet.addPoint(x, y, height);
    }
    const std::size_t bottomCenter = tet.addPoint(0.0, 0.0, 0.0);
    const std::size_t topCenter = tet.addPoint(0.0, 0.0, height);
    tet.tetrahedralize();
    AETHER_CHECK(tet.satisfiesDelaunayProperty());

    std::vector<std::array<std::size_t, 3>> facets;
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t j = (i + 1) % n;
        facets.push_back({bottomCenter, bottomRim[j], bottomRim[i]});
        facets.push_back({topCenter, topRim[i], topRim[j]});
        facets.push_back({bottomRim[i], bottomRim[j], topRim[j]});
        facets.push_back({bottomRim[i], topRim[j], topRim[i]});
    }
    AETHER_CHECK(facets.size() == 4 * n);
    AETHER_CHECK(!tet.missingFacets(facets).empty()); // measured: the ambiguity is real here too

    const auto result = tet.recoverFacets(facets, 4);
    AETHER_CHECK(result.recoveredFacets.size() == facets.size());
    AETHER_CHECK(result.unrecovered.empty()); // measured: the flip alone reaches every facet
    AETHER_CHECK(tet.satisfiesDelaunayProperty());

    const double baseArea = static_cast<double>(n) * 0.5 * radius * radius * std::sin(2.0 * kPi / static_cast<double>(n));
    AETHER_CHECK(nearlyEqual(sumTetrahedraVolume(tet), baseArea * height, 1e-9));
}

// Cross-validates tetrahedralize() (point location + adjacency, O(cavity
// size) per inserted point) against tetrahedralizeReference() (the
// original full-scan Bowyer-Watson, O(current tetrahedron count) per
// point) on the same points added in the same order. Delaunay
// tetrahedralization is unique for a point set once co-spherical ties are
// broken -- and both algorithms break them the same way, via the identical
// exact orientation3D/inSphere3D predicates -- so the two must produce
// *exactly* the same set of tetrahedra, not just the same count or the
// same total volume. Insertion order into each algorithm's own internal
// vector may legitimately differ, so tetrahedra are compared as a set of
// sorted 4-vertex-index tuples, not as an ordered sequence.
void checkFastMatchesReference(const std::vector<Vector3>& points, const char* label) {
    std::printf("  [mesh_tests] tetrahedralize() vs tetrahedralizeReference(): %s (%zu pontos)\n", label,
                points.size());

    DelaunayTetrahedralization3D fast;
    DelaunayTetrahedralization3D reference;
    for (const Vector3& p : points) {
        fast.addPoint(p.x, p.y, p.z);
        reference.addPoint(p.x, p.y, p.z);
    }
    fast.tetrahedralize();
    reference.tetrahedralizeReference();

    AETHER_CHECK(fast.satisfiesDelaunayProperty());
    AETHER_CHECK(reference.satisfiesDelaunayProperty());
    AETHER_CHECK(fast.tetrahedronCount() == reference.tetrahedronCount());
    AETHER_CHECK(nearlyEqual(sumTetrahedraVolume(fast), sumTetrahedraVolume(reference), 1e-9));

    std::set<std::array<std::size_t, 4>> fastSet;
    for (std::size_t t = 0; t < fast.tetrahedronCount(); ++t) {
        std::array<std::size_t, 4> v = fast.tetrahedron(t).vertices;
        std::sort(v.begin(), v.end());
        fastSet.insert(v);
    }
    std::set<std::array<std::size_t, 4>> referenceSet;
    for (std::size_t t = 0; t < reference.tetrahedronCount(); ++t) {
        std::array<std::size_t, 4> v = reference.tetrahedron(t).vertices;
        std::sort(v.begin(), v.end());
        referenceSet.insert(v);
    }
    AETHER_CHECK(fastSet == referenceSet);
}

// Reuses the exact point constructions of every geometrically tricky case
// this class already has a dedicated test for -- the hull-adjacent sliver
// that once broke on a too-small super-tetrahedron, the sparse cube, the
// dense prism -- plus two fresh jittered lattices for insertion-order/
// cavity-size diversity beyond the catalogued fixtures.
void testFastTetrahedralizeMatchesReferenceImplementation() {
    // Jittered 3x3x3 lattice (27 points), same formula as
    // testDelaunayTetrahedralizationSatisfiesDelaunayProperty().
    {
        std::vector<Vector3> points;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    const double x = static_cast<double>(i) + 0.1 * std::sin(3.7 * i + 1.3 * j + 2.4 * k);
                    const double y = static_cast<double>(j) + 0.1 * std::cos(2.1 * i + 4.2 * j + 1.1 * k);
                    const double z = static_cast<double>(k) + 0.1 * std::sin(1.6 * i + 2.9 * j + 3.3 * k);
                    points.emplace_back(x, y, z);
                }
            }
        }
        checkFastMatchesReference(points, "rede jitterada 3x3x3");
    }

    // Square bipyramid (6 points), same as
    // testDelaunayTetrahedralizationOfBipyramidPartitionsExactVolume().
    checkFastMatchesReference({Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0), Vector3(-1.0, 0.0, 0.0),
                                Vector3(0.0, -1.0, 0.0), Vector3(0.0, 0.0, 2.0), Vector3(0.0, 0.0, -1.0)},
                               "bipiramide quadrada");

    // The 7-point hull-adjacent sliver case, verbatim, same as
    // testTetrahedralizationFillsConvexHullOnHullAdjacentSliverCase().
    checkFastMatchesReference(
        {Vector3(1.4737617982207118, -0.8528567809098169, -0.09741757121024008),
         Vector3(0.8818641874977176, -1.5726403556472421, -1.5119823378293402),
         Vector3(1.7058871049329896, -0.5884653761446414, 0.4584386055893468),
         Vector3(-0.21987487313532084, -0.5879755232108659, -0.548828267652588),
         Vector3(-0.4376550633723242, -1.8110581604919749, 0.034080232074870764),
         Vector3(1.540759526282839, -1.6205730261557352, 0.6577405728359289),
         Vector3(0.9445458377233287, 1.12204139437512, -1.1761955929033237)},
        "sliver perto do casco (7 pontos)");

    // Sparse axis-aligned cube (8 points, no interior points), same as
    // testTetrahedralization3DFlipRecoversSomeCoplanarQuadFacetsHonestlyReportsTheRest().
    checkFastMatchesReference({Vector3(0.0, 0.0, 0.0), Vector3(2.0, 0.0, 0.0), Vector3(2.0, 2.0, 0.0),
                                Vector3(0.0, 2.0, 0.0), Vector3(0.0, 0.0, 2.0), Vector3(2.0, 0.0, 2.0),
                                Vector3(2.0, 2.0, 2.0), Vector3(0.0, 2.0, 2.0)},
                               "cubo esparso (8 pontos)");

    // Hollow double octahedron (12 points), same as
    // testTetrahedralization3DRecoversHollowOctahedronFacetsAndCarvesHoleWithExactVolume().
    {
        std::vector<Vector3> points;
        const Vector3 center(1.0, 1.0, 1.0);
        for (double radius : {2.0, 0.5}) {
            points.push_back(center + Vector3(radius, 0.0, 0.0));
            points.push_back(center + Vector3(-radius, 0.0, 0.0));
            points.push_back(center + Vector3(0.0, radius, 0.0));
            points.push_back(center + Vector3(0.0, -radius, 0.0));
            points.push_back(center + Vector3(0.0, 0.0, radius));
            points.push_back(center + Vector3(0.0, 0.0, -radius));
        }
        checkFastMatchesReference(points, "octaedro oco duplo (12 pontos)");
    }

    // Dense 24-sided polygonal prism (50 points), same as
    // testTetrahedralization3DRecoversPolygonalPrismSideFacetsLikeACylinder().
    {
        constexpr std::size_t n = 24;
        constexpr double radius = 0.7;
        constexpr double height = 5.3;
        constexpr double kPi = 3.14159265358979323846;
        std::vector<Vector3> points;
        for (std::size_t i = 0; i < n; ++i) {
            const double angle = 2.0 * kPi * static_cast<double>(i) / static_cast<double>(n);
            const double x = radius * std::cos(angle);
            const double y = radius * std::sin(angle);
            points.emplace_back(x, y, 0.0);
            points.emplace_back(x, y, height);
        }
        points.emplace_back(0.0, 0.0, 0.0);
        points.emplace_back(0.0, 0.0, height);
        checkFastMatchesReference(points, "prisma poligonal denso (50 pontos)");
    }

    // Two fresh, larger jittered lattices, for insertion-order and
    // cavity-size diversity beyond the fixtures above -- same jitter
    // formula, larger n.
    for (std::size_t n : {4, 6}) {
        std::vector<Vector3> points;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t k = 0; k < n; ++k) {
                    const double fi = static_cast<double>(i);
                    const double fj = static_cast<double>(j);
                    const double fk = static_cast<double>(k);
                    const double x = fi + 0.1 * std::sin(3.7 * fi + 1.3 * fj + 2.4 * fk);
                    const double y = fj + 0.1 * std::cos(2.1 * fi + 4.2 * fj + 1.1 * fk);
                    const double z = fk + 0.1 * std::sin(1.6 * fi + 2.9 * fj + 3.3 * fk);
                    points.emplace_back(x, y, z);
                }
            }
        }
        char label[64];
        std::snprintf(label, sizeof(label), "rede jitterada %zux%zux%zu", n, n, n);
        checkFastMatchesReference(points, label);
    }
}

// ---------------------------------------------------------------------------
// TetrahedralMesh (Fase 2.1 do ROADMAP): face connectivity, the prerequisite
// for any unstructured finite-volume solver.
//
// Every check below is an *identity*, not a tolerance on physics: these hold
// for any valid tetrahedral mesh regardless of what is later solved on it,
// which is what makes them able to catch the two bugs this code could
// plausibly have (a misoriented face normal, or a face registered to only one
// of the two cells that share it).
// ---------------------------------------------------------------------------

// A deterministic, genuinely 3D point cloud: a unit cube's 8 corners plus
// interior points on a jittered lattice. Jittered rather than regular
// because a perfect lattice is massively co-spherical, which is Delaunay's
// degenerate tie case -- a worse test than a validation, exactly as
// testDelaunayTriangulation2D already notes for the 2D case.
DelaunayTetrahedralization3D buildJitteredCubeTetrahedralization(std::size_t interiorPerAxis) {
    DelaunayTetrahedralization3D tets;
    for (int corner = 0; corner < 8; ++corner) {
        tets.addPoint((corner & 1) ? 1.0 : 0.0, (corner & 2) ? 1.0 : 0.0, (corner & 4) ? 1.0 : 0.0);
    }
    for (std::size_t i = 1; i <= interiorPerAxis; ++i) {
        for (std::size_t j = 1; j <= interiorPerAxis; ++j) {
            for (std::size_t k = 1; k <= interiorPerAxis; ++k) {
                const double step = 1.0 / static_cast<double>(interiorPerAxis + 1);
                const double x = static_cast<double>(i) * step;
                const double y = static_cast<double>(j) * step;
                const double z = static_cast<double>(k) * step;
                const double jitter = 0.13 * step;
                tets.addPoint(x + jitter * std::sin(7.0 * x + 3.0 * y),
                              y + jitter * std::sin(5.0 * y + 2.0 * z),
                              z + jitter * std::sin(3.0 * z + 4.0 * x));
            }
        }
    }
    tets.tetrahedralize();
    return tets;
}

// **The discrete divergence theorem.** For any closed polyhedron the sum of
// its outward face area vectors is identically zero. This is the single
// strongest check available on face connectivity: it fails immediately if
// any face normal points the wrong way, if a face is missing from a cell's
// list, or if the owner/neighbour orientation convention is inconsistent --
// and it holds independently of mesh quality, so it needs no tolerance
// beyond roundoff.
void testTetrahedralMeshClosesEveryCell() {
    const DelaunayTetrahedralization3D tets = buildJitteredCubeTetrahedralization(3);
    const TetrahedralMesh mesh = TetrahedralMesh::fromTetrahedralization(tets);

    AETHER_CHECK(mesh.cellCount() > 0);

    double worst = 0.0;
    double smallestCellScale = 1e300;
    for (std::size_t cell = 0; cell < mesh.cellCount(); ++cell) {
        AETHER_CHECK(mesh.cellFaces(cell).size() == 4); // a tetrahedron has exactly four faces
        worst = std::max(worst, mesh.cellAreaVectorSum(cell).norm());
        // Normalizing by a cell's own area scale is what makes "zero"
        // meaningful: an absolute bound would be trivially passed by a tiny
        // sliver and unfairly strict on a large cell.
        double areaScale = 0.0;
        for (std::size_t faceIndex : mesh.cellFaces(cell)) {
            areaScale += mesh.outwardAreaVector(cell, faceIndex).norm();
        }
        smallestCellScale = std::min(smallestCellScale, areaScale);
    }
    std::printf("  [mesh_tests] TetrahedralMesh: %zu celulas, pior |soma dos vetores de area| = %.3e\n",
                mesh.cellCount(), worst);
    AETHER_CHECK(worst < 1e-12 * smallestCellScale + 1e-14);
}

// Euler-style bookkeeping: every tetrahedron contributes four faces, and
// each interior face is contributed by exactly two of them while each
// boundary face is contributed by one. So 4*cells = 2*interior + boundary
// must hold exactly, in integers. This catches a face that was created
// twice (hash collision or key-ordering bug) or never matched up.
void testTetrahedralMeshFaceCountsBalance() {
    const DelaunayTetrahedralization3D tets = buildJitteredCubeTetrahedralization(3);
    const TetrahedralMesh mesh = TetrahedralMesh::fromTetrahedralization(tets);

    const std::size_t boundary = mesh.boundaryFaceCount();
    const std::size_t interior = mesh.faceCount() - boundary;
    std::printf("  [mesh_tests] TetrahedralMesh: %zu faces (%zu internas, %zu de contorno)\n",
                mesh.faceCount(), interior, boundary);

    AETHER_CHECK(4 * mesh.cellCount() == 2 * interior + boundary);
    AETHER_CHECK(boundary > 0); // a bounded domain must have a boundary
    AETHER_CHECK(interior > 0); // and this mesh is far too refined to have none

    // Owner/neighbour must be a genuine pairing, never self-referential.
    for (std::size_t f = 0; f < mesh.faceCount(); ++f) {
        const auto& face = mesh.face(f);
        AETHER_CHECK(face.owner < mesh.cellCount());
        if (!mesh.isBoundaryFace(f)) {
            AETHER_CHECK(face.neighbour < mesh.cellCount());
            AETHER_CHECK(face.neighbour != face.owner);
        }
    }
}

// The tetrahedralization of a point set fills exactly its convex hull, and
// here that hull is the unit cube (its eight corners are among the input
// points, and every other point is inside). So the cell volumes must sum to
// exactly 1 -- an independent check of the volume formula and of the claim
// that the cells tile the domain without gaps or overlaps.
void testTetrahedralMeshFillsTheCubeExactly() {
    const DelaunayTetrahedralization3D tets = buildJitteredCubeTetrahedralization(3);
    const TetrahedralMesh mesh = TetrahedralMesh::fromTetrahedralization(tets);

    const double volume = mesh.totalVolume();
    std::printf("  [mesh_tests] TetrahedralMesh: volume total = %.15f (esperado 1.0)\n", volume);
    std::fflush(stdout);
    AETHER_CHECK(nearlyEqual(volume, 1.0, 1e-12));
}

// The two sides of an interior face must see exactly opposite area vectors,
// or flux would not be conserved between neighbours -- the property an FVM
// discretization relies on to be conservative at all. Checked as exact
// negation, since outwardAreaVector() flips one stored vector rather than
// recomputing it.
void testTetrahedralMeshFaceOrientationIsAntisymmetric() {
    const DelaunayTetrahedralization3D tets = buildJitteredCubeTetrahedralization(2);
    const TetrahedralMesh mesh = TetrahedralMesh::fromTetrahedralization(tets);

    std::size_t checked = 0;
    for (std::size_t f = 0; f < mesh.faceCount(); ++f) {
        if (mesh.isBoundaryFace(f)) {
            continue;
        }
        const auto& face = mesh.face(f);
        const Vector3 fromOwner = mesh.outwardAreaVector(face.owner, f);
        const Vector3 fromNeighbour = mesh.outwardAreaVector(face.neighbour, f);
        AETHER_CHECK(fromOwner.x == -fromNeighbour.x);
        AETHER_CHECK(fromOwner.y == -fromNeighbour.y);
        AETHER_CHECK(fromOwner.z == -fromNeighbour.z);
        ++checked;
    }
    AETHER_CHECK(checked > 0);
}

// The area vector must point from the owner towards the neighbour, not
// merely along that line. Getting this backwards on some faces would still
// satisfy the closure identity above (it sums over a cell, and a flipped
// pair cancels within it), so it needs its own check -- otherwise every
// flux computed later would carry the wrong sign on those faces.
void testTetrahedralMeshAreaVectorsPointOwnerToNeighbour() {
    const DelaunayTetrahedralization3D tets = buildJitteredCubeTetrahedralization(2);
    const TetrahedralMesh mesh = TetrahedralMesh::fromTetrahedralization(tets);

    std::size_t checked = 0;
    for (std::size_t f = 0; f < mesh.faceCount(); ++f) {
        if (mesh.isBoundaryFace(f)) {
            continue;
        }
        const auto& face = mesh.face(f);
        const Vector3 ownerToNeighbour = mesh.cellCentroid(face.neighbour) - mesh.cellCentroid(face.owner);
        AETHER_CHECK(face.areaVector.dot(ownerToNeighbour) > 0.0);
        ++checked;
    }
    std::printf("  [mesh_tests] TetrahedralMesh: %zu faces internas com orientacao owner->neighbour ok\n",
                checked);
    std::fflush(stdout);
    AETHER_CHECK(checked > 0);
}

// **The claim the floating-point filter makes is that it changes no
// answer** -- it only decides whether the exact BigInt path can be
// skipped. That is exactly the kind of claim this project refuses to take
// on faith, so it is checked directly: run both paths over the same
// inputs and require the signs to be identical every time.
//
// The inputs are chosen to include the cases where a filter could
// plausibly be wrong, not just easy ones. Random points exercise the fast
// path (and confirm it really is taken, via the fallback counts printed
// below -- a filter that never fires would pass this test while being
// useless). The degenerate families are the adversarial half: four points
// forced exactly coplanar, and a fifth point placed exactly on the
// circumsphere of the other four, are precisely the configurations where
// the true determinant is zero and any error bound that is even slightly
// too tight would confidently return a nonzero sign.
void testPredicateFilterAgreesWithExactArithmetic() {
    std::uint64_t state = 0x9E3779B97F4A7C15ull;
    const auto nextDouble = [&state]() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return static_cast<double>(state >> 11) / 9007199254740992.0; // [0,1)
    };
    const auto randomPoint = [&nextDouble]() {
        return Vector3{nextDouble() * 2.0 - 1.0, nextDouble() * 2.0 - 1.0, nextDouble() * 2.0 - 1.0};
    };

    std::size_t orientationChecks = 0;
    std::size_t orientationFallbacks = 0;
    std::size_t inSphereChecks = 0;
    std::size_t inSphereFallbacks = 0;

    const auto checkOrientation = [&](const Vector3& a, const Vector3& b, const Vector3& c,
                                       const Vector3& d) {
        const int exact = orientation3DExact(a, b, c, d);
        AETHER_CHECK(orientation3D(a, b, c, d) == exact);
        ++orientationChecks;
        if (exact == 0) {
            ++orientationFallbacks; // a zero sign can only come from the exact path
        }
    };
    const auto checkInSphere = [&](const Vector3& a, const Vector3& b, const Vector3& c,
                                    const Vector3& d, const Vector3& p) {
        const int exact = inSphere3DExact(a, b, c, d, p);
        AETHER_CHECK(inSphere3D(a, b, c, d, p) == exact);
        ++inSphereChecks;
        if (exact == 0) {
            ++inSphereFallbacks;
        }
    };

    for (int i = 0; i < 20000; ++i) {
        const Vector3 a = randomPoint();
        const Vector3 b = randomPoint();
        const Vector3 c = randomPoint();
        const Vector3 d = randomPoint();
        const Vector3 p = randomPoint();
        checkOrientation(a, b, c, d);
        checkInSphere(a, b, c, d, p);

        // Exactly coplanar: d is built from a, b, c by an exact affine
        // combination with power-of-two weights, so no rounding creeps in
        // and the true orientation determinant is exactly zero.
        const Vector3 coplanar{a.x + 0.5 * (b.x - a.x) + 0.25 * (c.x - a.x),
                                a.y + 0.5 * (b.y - a.y) + 0.25 * (c.y - a.y),
                                a.z + 0.5 * (b.z - a.z) + 0.25 * (c.z - a.z)};
        checkOrientation(a, b, c, coplanar);

        // Exactly co-spherical: five points on a sphere of radius 1 about
        // the origin, read off axis directions and their exact negations,
        // so the in-sphere determinant is exactly zero by construction.
        const Vector3 s0{1.0, 0.0, 0.0};
        const Vector3 s1{-1.0, 0.0, 0.0};
        const Vector3 s2{0.0, 1.0, 0.0};
        const Vector3 s3{0.0, 0.0, 1.0};
        const Vector3 s4{0.0, -1.0, 0.0};
        checkInSphere(s0, s2, s3, s1, s4);

        // A sliver: three nearly-collinear points plus a fourth, the
        // near-degenerate shape the class header records as having caused
        // a real bug in the unfiltered double era.
        const double squash = 1e-11;
        const Vector3 t0{0.0, 0.0, 0.0};
        const Vector3 t1{1.0, squash * nextDouble(), 0.0};
        const Vector3 t2{2.0, -squash * nextDouble(), 0.0};
        const Vector3 t3{1.0, 0.0, 1.0};
        checkOrientation(t0, t1, t2, t3);
        checkInSphere(t0, t1, t2, t3, randomPoint());
    }

    std::printf("  [mesh_tests] filtro de predicados: %zu orientacoes e %zu in-sphere conferem com "
                "aritmetica exata (%zu e %zu degenerados de sinal zero)\n",
                orientationChecks, inSphereChecks, orientationFallbacks, inSphereFallbacks);
    std::fflush(stdout);
    // Degenerate cases must actually have been produced -- otherwise this
    // test would be measuring only the easy half of the input space.
    AETHER_CHECK(orientationFallbacks > 0);
    AETHER_CHECK(inSphereFallbacks > 0);
}

} // namespace

int main() {
    testTetrahedralMeshClosesEveryCell();
    testTetrahedralMeshFaceCountsBalance();
    testTetrahedralMeshFillsTheCubeExactly();
    testTetrahedralMeshFaceOrientationIsAntisymmetric();
    testTetrahedralMeshAreaVectorsPointOwnerToNeighbour();
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
    testTetrahedralization3DFlipRecoversSomeCoplanarQuadFacetsHonestlyReportsTheRest();
    testTetrahedralization3DRecoversPolygonalPrismSideFacetsLikeACylinder();
    testFastTetrahedralizeMatchesReferenceImplementation();
    testPredicateFilterAgreesWithExactArithmetic();
    std::puts("aether_mesh_tests: all tests passed");
    return 0;
}
