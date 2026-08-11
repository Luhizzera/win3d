#include "aether/geometry/ObjIO.hpp"
#include "aether/geometry/StlIO.hpp"
#include "aether/geometry/TriangleMesh.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>

using namespace aether::core;
using namespace aether::geometry;

namespace {

bool nearlyEqual(double a, double b, double tolerance = 1e-6) {
    return std::fabs(a - b) <= tolerance;
}

// Right tetrahedron with one vertex at the origin and the others at the
// unit points on each axis. Volume = 1/6, easy to check by hand.
TriangleMesh buildUnitTetrahedron() {
    TriangleMesh mesh;
    const std::size_t v0 = mesh.addVertex({0.0, 0.0, 0.0});
    const std::size_t v1 = mesh.addVertex({1.0, 0.0, 0.0});
    const std::size_t v2 = mesh.addVertex({0.0, 1.0, 0.0});
    const std::size_t v3 = mesh.addVertex({0.0, 0.0, 1.0});

    mesh.addTriangle(v0, v2, v1);
    mesh.addTriangle(v0, v1, v3);
    mesh.addTriangle(v0, v3, v2);
    mesh.addTriangle(v1, v2, v3);
    return mesh;
}

void testVolumeAndArea() {
    const TriangleMesh mesh = buildUnitTetrahedron();
    AETHER_CHECK(nearlyEqual(mesh.volume(), 1.0 / 6.0));

    const double expectedArea = 1.5 + std::sqrt(3.0) / 2.0;
    AETHER_CHECK(nearlyEqual(mesh.surfaceArea(), expectedArea));
}

void testWatertightness() {
    const TriangleMesh closed = buildUnitTetrahedron();
    AETHER_CHECK(closed.isWatertight());
    AETHER_CHECK(closed.findBoundaryEdges().empty());

    TriangleMesh open;
    const std::size_t v0 = open.addVertex({0.0, 0.0, 0.0});
    const std::size_t v1 = open.addVertex({1.0, 0.0, 0.0});
    const std::size_t v2 = open.addVertex({0.0, 1.0, 0.0});
    const std::size_t v3 = open.addVertex({0.0, 0.0, 1.0});
    open.addTriangle(v0, v2, v1);
    open.addTriangle(v0, v1, v3);
    open.addTriangle(v0, v3, v2);
    // Deliberately omit the (v1, v2, v3) face to leave a hole.

    AETHER_CHECK(!open.isWatertight());
    AETHER_CHECK(open.findBoundaryEdges().size() == 3);
}

void testWeldVertices() {
    // Two triangles sharing an edge, stored the way STL would: every vertex
    // duplicated per triangle, no shared indices.
    TriangleMesh mesh;
    const std::size_t a0 = mesh.addVertex({0.0, 0.0, 0.0});
    const std::size_t a1 = mesh.addVertex({1.0, 0.0, 0.0});
    const std::size_t a2 = mesh.addVertex({0.0, 1.0, 0.0});
    mesh.addTriangle(a0, a1, a2);

    const std::size_t b0 = mesh.addVertex({1.0, 0.0, 0.0});
    const std::size_t b1 = mesh.addVertex({1.0, 1.0, 0.0});
    const std::size_t b2 = mesh.addVertex({0.0, 1.0, 0.0});
    mesh.addTriangle(b0, b1, b2);

    AETHER_CHECK(mesh.vertexCount() == 6);
    const std::size_t merged = mesh.weldVertices();
    AETHER_CHECK(merged == 2);
    AETHER_CHECK(mesh.vertexCount() == 4);
    AETHER_CHECK(!mesh.isWatertight()); // open quad: 4 boundary edges
    AETHER_CHECK(mesh.findBoundaryEdges().size() == 4);
}

void testRemoveDegenerateTriangles() {
    TriangleMesh mesh = buildUnitTetrahedron();
    const std::size_t v0 = mesh.addVertex({5.0, 5.0, 5.0});
    mesh.addTriangle(v0, v0, v0); // fully degenerate: repeated vertex

    AETHER_CHECK(mesh.triangleCount() == 5);
    const std::size_t removed = mesh.removeDegenerateTriangles();
    AETHER_CHECK(removed == 1);
    AETHER_CHECK(mesh.triangleCount() == 4);
}

void testStlRoundTrip() {
    const TriangleMesh original = buildUnitTetrahedron();
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_geometry_test_tetra.stl";

    saveStlBinary(original, path.string());
    const TriangleMesh loaded = loadStl(path.string());

    AETHER_CHECK(loaded.vertexCount() == 4);
    AETHER_CHECK(loaded.triangleCount() == 4);
    AETHER_CHECK(nearlyEqual(loaded.volume(), 1.0 / 6.0));
    AETHER_CHECK(loaded.isWatertight());

    std::filesystem::remove(path);
}

void testObjRoundTrip() {
    const TriangleMesh original = buildUnitTetrahedron();
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_geometry_test_tetra.obj";

    saveObj(original, path.string());
    const TriangleMesh loaded = loadObj(path.string());

    AETHER_CHECK(loaded.vertexCount() == 4);
    AETHER_CHECK(loaded.triangleCount() == 4);
    AETHER_CHECK(nearlyEqual(loaded.volume(), 1.0 / 6.0));
    AETHER_CHECK(loaded.isWatertight());

    std::filesystem::remove(path);
}

// A hand-written OBJ face with more than 3 vertices (a planar quad, the
// common case real exporters produce) must be fan-triangulated into 2
// triangles from its first vertex, not silently dropped or mis-parsed.
// Also exercises the "v/vt/vn" face-token form (only the leading vertex
// index should be used) since that is what most real-world OBJ exporters
// emit, not the bare "v"-only form saveObj() itself writes.
void testObjQuadFaceIsFanTriangulated() {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_geometry_test_quad.obj";
    {
        std::ofstream file(path);
        file << "v 0.0 0.0 0.0\n";
        file << "v 2.0 0.0 0.0\n";
        file << "v 2.0 1.0 0.0\n";
        file << "v 0.0 1.0 0.0\n";
        file << "vn 0.0 0.0 1.0\n";
        file << "f 1/1/1 2/1/1 3/1/1 4/1/1\n";
    }

    const TriangleMesh mesh = loadObj(path.string());
    AETHER_CHECK(mesh.vertexCount() == 4);
    AETHER_CHECK(mesh.triangleCount() == 2); // fan-triangulated quad: n-2 triangles

    double summedArea = 0.0;
    for (std::size_t t = 0; t < mesh.triangleCount(); ++t) {
        summedArea += mesh.faceArea(t);
    }
    AETHER_CHECK(nearlyEqual(summedArea, 2.0)); // 2x1 rectangle

    std::filesystem::remove(path);
}

} // namespace

int main() {
    testVolumeAndArea();
    testWatertightness();
    testWeldVertices();
    testRemoveDegenerateTriangles();
    testStlRoundTrip();
    testObjRoundTrip();
    testObjQuadFaceIsFanTriangulated();
    std::puts("aether_geometry_tests: all tests passed");
    return 0;
}
