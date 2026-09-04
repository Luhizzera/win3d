#include "aether/geometry/ObjIO.hpp"
#include "aether/geometry/StepIO.hpp"
#include "aether/geometry/StlIO.hpp"
#include "aether/geometry/TriangleMesh.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>

#ifdef AETHER_HAVE_OPENCASCADE
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <STEPControl_StepModelType.hxx>
#include <STEPControl_Writer.hxx>
#include <TopoDS_Shape.hxx>
#endif

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

// **Both tests below exercise the tokenizer's own leniency, not a
// property of STEP itself.** They hand-write only the entities StepIO.hpp
// says were checked against the schema -- no PRODUCT/PRODUCT_DEFINITION/
// SHAPE_DEFINITION_REPRESENTATION/GEOMETRIC_REPRESENTATION_CONTEXT
// boilerplate a real CAD export always carries -- because the tokenizer
// looks for a FACETED_BREP by type name directly and does not need any of
// it. A real CAD kernel does: measured directly, OpenCASCADE's
// STEPControl_Reader transfers *no shape at all* from this same minimal
// tetrahedron file ("file transferred no shape"), because it has no root
// entity of the kind its transfer mechanism looks for. That is not a
// defect in either implementation -- the tokenizer's whole point is
// accepting exactly the restricted grammar it documents, and a real
// kernel's insistence on a complete, standards-shaped file is what makes
// it trustworthy on real-world input. So both tests compile only in the
// no-OpenCASCADE configuration; testStepLoadsCurvedCylinderThroughOpenCascade
// below is the equivalent gate for the OpenCASCADE path, using a file
// OpenCASCADE's own writer produced rather than one hand-authored to the
// tokenizer's minimal grammar.
#ifndef AETHER_HAVE_OPENCASCADE

// A hand-written, minimal ISO-10303-21 file describing the same unit
// tetrahedron as buildUnitTetrahedron() above, as a FACETED_BREP: four
// CARTESIAN_POINTs, four POLY_LOOPs, four FACE_OUTER_BOUNDs, four FACEs, one
// CLOSED_SHELL. Every entity shape used here is exactly what StepIO.hpp's
// header comment says was checked against ISO 10303-42's published EXPRESS
// schema, so this file is deliberately not copied from any real CAD
// export -- it is the grammar's own minimal instantiation.
//
// One face (#21, the (v0,v1,v3) face) is deliberately written with its loop
// in the *wrong* rotational sense and orientation=.F., so the round trip
// also exercises FACE_BOUND's orientation flag actually reversing the loop,
// not just the common orientation=.T. case the other three faces use.
void testStepLoadsFacetedTetrahedron() {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_geometry_test_tetra.step";
    {
        std::ofstream file(path);
        file << "ISO-10303-21;\n";
        file << "HEADER;\n";
        file << "FILE_DESCRIPTION((''),'2;1');\n";
        file << "FILE_NAME('tetra.step','2024-01-01T00:00:00',(''),(''),'','','');\n";
        file << "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n";
        file << "ENDSEC;\n";
        file << "DATA;\n";
        file << "#1=CARTESIAN_POINT('',(0.0,0.0,0.0));\n";  // v0
        file << "#2=CARTESIAN_POINT('',(1.0,0.0,0.0));\n";  // v1
        file << "#3=CARTESIAN_POINT('',(0.0,1.0,0.0));\n";  // v2
        file << "#4=CARTESIAN_POINT('',(0.0,0.0,1.0));\n";  // v3
        file << "#10=POLY_LOOP('',(#1,#3,#2));\n";          // (v0,v2,v1), already outward
        file << "#11=POLY_LOOP('',(#1,#4,#2));\n";          // (v0,v3,v1): wrong sense on purpose
        file << "#12=POLY_LOOP('',(#1,#4,#3));\n";          // (v0,v3,v2), already outward
        file << "#13=POLY_LOOP('',(#2,#3,#4));\n";          // (v1,v2,v3), already outward
        file << "#20=FACE_OUTER_BOUND('',#10,.T.);\n";
        file << "#21=FACE_OUTER_BOUND('',#11,.F.);\n"; // reversal must fix #11's wrong sense
        file << "#22=FACE_OUTER_BOUND('',#12,.T.);\n";
        file << "#23=FACE_OUTER_BOUND('',#13,.T.);\n";
        file << "#30=FACE('',(#20));\n";
        file << "#31=FACE('',(#21));\n";
        file << "#32=FACE('',(#22));\n";
        file << "#33=FACE('',(#23));\n";
        file << "#40=CLOSED_SHELL('',(#30,#31,#32,#33));\n";
        file << "#50=FACETED_BREP('',#40);\n";
        file << "ENDSEC;\n";
        file << "END-ISO-10303-21;\n";
    }

    const StepLoadResult loaded = loadStep(path.string());
    AETHER_CHECK(loaded.unsupportedFeatures.empty());
    AETHER_CHECK(loaded.mesh.vertexCount() == 4);
    AETHER_CHECK(loaded.mesh.triangleCount() == 4);
    AETHER_CHECK(nearlyEqual(loaded.mesh.volume(), 1.0 / 6.0));

    // Every face must wind outward -- including #31, which only ends up
    // outward at all because the loader honoured orientation=.F. Checked
    // against the mesh's own centroid rather than hardcoded normals, so
    // this genuinely fails if the reversal is skipped or inverted.
    Vector3 centroid;
    for (std::size_t v = 0; v < loaded.mesh.vertexCount(); ++v) {
        centroid += loaded.mesh.vertex(v);
    }
    centroid = centroid / static_cast<double>(loaded.mesh.vertexCount());

    for (std::size_t t = 0; t < loaded.mesh.triangleCount(); ++t) {
        const Triangle& tri = loaded.mesh.triangle(t);
        Vector3 faceCentroid;
        for (const std::size_t vertexIndex : tri.vertices) {
            faceCentroid += loaded.mesh.vertex(vertexIndex);
        }
        faceCentroid = faceCentroid / 3.0;
        AETHER_CHECK(loaded.mesh.faceNormal(t).dot(faceCentroid - centroid) > 0.0);
    }

    std::filesystem::remove(path);
}

// A FACE whose bound is an EDGE_LOOP rather than a POLY_LOOP is exactly the
// shape a curved (non-faceted) boundary takes in a real STEP file. The
// loader must name it in unsupportedFeatures rather than silently produce
// an empty or wrong mesh -- StepLoadResult's whole reason to carry that
// field rather than just throwing or returning nothing.
// Tests the *tokenizer's own* unsupported-feature-reporting contract
// (StepIO.cpp's no-OpenCASCADE path): a structurally-parseable-but-
// dangling-reference face is recognized as not faceted and reported by
// name rather than crashing or silently producing an empty/wrong mesh.
// Compiled only in this build configuration because the function under
// test -- the tokenizer -- is not the one loadStep() calls once
// AETHER_HAVE_OPENCASCADE is defined (a real CAD kernel reads the file
// then, and has no equivalent "valid syntax, unsupported geometry"
// category to report the same way -- see StepIO.hpp's own comment).
void testStepReportsNonFacetedBoundaryAsUnsupported() {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_geometry_test_curved.step";
    {
        std::ofstream file(path);
        file << "ISO-10303-21;\n";
        file << "HEADER;\n";
        file << "FILE_DESCRIPTION((''),'2;1');\n";
        file << "FILE_NAME('curved.step','2024-01-01T00:00:00',(''),(''),'','','');\n";
        file << "FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\n";
        file << "ENDSEC;\n";
        file << "DATA;\n";
        // #101..#103 (edges) are deliberately never defined: resolveLoop()
        // rejects this loop by its EDGE_LOOP type name before ever
        // dereferencing them, so a realistic EDGE_CURVE/VERTEX_POINT chain
        // is not needed just to exercise that rejection.
        file << "#10=EDGE_LOOP('',(#101,#102,#103));\n";
        file << "#20=FACE_OUTER_BOUND('',#10,.T.);\n";
        file << "#30=FACE('',(#20));\n";
        file << "#40=CLOSED_SHELL('',(#30));\n";
        file << "#50=FACETED_BREP('',#40);\n";
        file << "ENDSEC;\n";
        file << "END-ISO-10303-21;\n";
    }

    const StepLoadResult loaded = loadStep(path.string());
    AETHER_CHECK(loaded.mesh.triangleCount() == 0);
    AETHER_CHECK(!loaded.unsupportedFeatures.empty());
    bool mentionsFace30 = false;
    for (const std::string& feature : loaded.unsupportedFeatures) {
        if (feature.find("#30") != std::string::npos) {
            mentionsFace30 = true;
        }
    }
    AETHER_CHECK(mentionsFace30);

    std::filesystem::remove(path);
}
#endif // !AETHER_HAVE_OPENCASCADE

#ifdef AETHER_HAVE_OPENCASCADE
// **The gate for ROADMAP Fase 5's curved-surface case.** A hand-authored
// Part 21 file cannot practically express a curved face the way
// testStepLoadsFacetedTetrahedron does for a planar one, so this builds
// the reference solid *with OpenCASCADE itself* -- a real cylinder,
// genuinely curved -- writes it to a real STEP file, and reads it back
// through loadStep()'s own public API exactly as any caller would: a
// real round trip through curved-surface STEP, still self-contained (no
// external fixture file).
void testStepLoadsCurvedCylinderThroughOpenCascade() {
    constexpr double kPi = 3.14159265358979323846;
    constexpr double kRadius = 2.0;
    constexpr double kHeight = 5.0;
    const TopoDS_Shape cylinder = BRepPrimAPI_MakeCylinder(kRadius, kHeight).Shape();

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_geometry_test_cylinder.step";
    {
        STEPControl_Writer writer;
        writer.Transfer(cylinder, STEPControl_AsIs);
        writer.Write(path.string().c_str());
    }

    const StepLoadResult loaded = loadStep(path.string());
    std::filesystem::remove(path);

    AETHER_CHECK(loaded.unsupportedFeatures.empty());
    AETHER_CHECK(loaded.mesh.triangleCount() > 0);

    const double exactVolume = kPi * kRadius * kRadius * kHeight;
    const double measuredVolume = loaded.mesh.volume();
    const double relativeError = std::fabs(measuredVolume - exactVolume) / exactVolume;
    std::printf("  [geometry_tests] cilindro curvo via OpenCASCADE: volume exato=%.6f medido=%.6f "
                "erro relativo=%.6e (%zu triangulos)\n",
                exactVolume, measuredVolume, relativeError, loaded.mesh.triangleCount());
    std::fflush(stdout);

    // Tessellation of a curved surface is inherently approximate.
    // Measured 8.30e-04 against the deflection tessellateShape()
    // (StepIO.cpp) actually uses (1e-3 relative linear, 0.3 rad angular) --
    // not assumed tight, the same distinction the Sod shock tube and
    // composite-wall tests draw elsewhere in this project. Set with
    // headroom above the measured value.
    AETHER_CHECK(relativeError < 2e-3);
}
#endif // AETHER_HAVE_OPENCASCADE

} // namespace

int main() {
    testVolumeAndArea();
    testWatertightness();
    testWeldVertices();
    testRemoveDegenerateTriangles();
    testStlRoundTrip();
    testObjRoundTrip();
    testObjQuadFaceIsFanTriangulated();
#ifndef AETHER_HAVE_OPENCASCADE
    testStepLoadsFacetedTetrahedron();
    testStepReportsNonFacetedBoundaryAsUnsupported();
#endif
#ifdef AETHER_HAVE_OPENCASCADE
    testStepLoadsCurvedCylinderThroughOpenCascade();
#endif
    std::puts("aether_geometry_tests: all tests passed");
    return 0;
}
