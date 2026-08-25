#include "aether/postprocessing/MarchingCubes3D.hpp"
#include "aether/postprocessing/MarchingSquares2D.hpp"
#include "aether/postprocessing/Streamline2D.hpp"
#include "aether/postprocessing/VtkWriter.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace aether::core;
using namespace aether::mesh;
using namespace aether::postprocessing;

namespace {

bool nearlyEqual(double a, double b, double tolerance) { return std::fabs(a - b) <= tolerance; }

// Streamline integration validated against the Taylor-Green vortex's
// stream function: for u = U0*cos(x)*sin(y), v = -U0*sin(x)*cos(y) (the
// same field TaylorGreenVortexSolver2D's own test validates as an exact
// solution of Navier-Stokes), psi = -U0*cos(x)*cos(y) is the exact stream
// function (du/dy=psi, -dv/dx=... standard check: dpsi/dy=u, -dpsi/dx=v,
// confirmed by direct differentiation). A streamline is, by definition, a
// curve of constant psi -- so integrating one via Streamline2D and
// checking psi stays constant along the whole path is a strong, self-
// contained, exact invariant check (the same style used throughout this
// project: an analytically-derivable property, not a recalled benchmark).
void testStreamlineConservesTaylorGreenStreamFunction() {
    const double kPi = 3.14159265358979323846;
    const std::size_t n = 64;
    const double length = 2.0 * kPi;
    const double u0 = 1.0;
    const double dx = length / static_cast<double>(n);

    std::vector<double> u(n * n);
    std::vector<double> v(n * n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * dx;
            const double y = (static_cast<double>(j) + 0.5) * dx;
            u[i + j * n] = u0 * std::cos(x) * std::sin(y);
            v[i + j * n] = -u0 * std::sin(x) * std::cos(y);
        }
    }

    Streamline2D tracer(n, n, length, length, u, v, /*periodic=*/true);

    const double x0 = 1.0;
    const double y0 = 0.5;
    const double psi0 = -u0 * std::cos(x0) * std::cos(y0);

    const auto path = tracer.trace(x0, y0, 0.01, 3000);
    AETHER_CHECK(path.size() > 1000); // did not immediately stagnate or exit

    double maxPsiDrift = 0.0;
    for (const Vector3& p : path) {
        const double psi = -u0 * std::cos(p.x) * std::cos(p.y);
        maxPsiDrift = std::max(maxPsiDrift, std::fabs(psi - psi0));
    }
    // Measured directly (not guessed): RK4 drift over 3000 steps at this
    // resolution/step size is ~7.7e-7; tolerance below keeps a comfortable
    // margin above that.
    AETHER_CHECK(maxPsiDrift < 1e-4);
}

// Marching squares validated against a field with a known exact contour:
// f(x,y) = x^2 + y^2 has the circle of radius r as its f = r^2 iso-
// contour, for any r. Every extracted segment endpoint must lie within
// ordinary grid-resolution error of that circle, and the total contour
// length must approximate its circumference (2*pi*r) -- both self-
// derivable from the field definition, not a recalled reference shape.
void testMarchingSquaresExtractsCircle() {
    const std::size_t n = 128;
    const double domain = 4.0; // covers [-2, 2] in both x and y
    const double dx = domain / static_cast<double>(n);
    const double radius = 1.0;
    const double isoValue = radius * radius;

    std::vector<double> field(n * n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * dx - domain / 2.0;
            const double y = (static_cast<double>(j) + 0.5) * dx - domain / 2.0;
            field[i + j * n] = x * x + y * y;
        }
    }

    const auto segments = marchingSquares2D(n, n, domain, domain, field, isoValue);
    AETHER_CHECK(!segments.empty());

    double maxRadiusError = 0.0;
    double totalLength = 0.0;
    for (const Segment2D& seg : segments) {
        // Endpoints are on the (shifted) grid; recenter before checking
        // distance from the circle's actual center.
        const double ax = seg.a.x - domain / 2.0;
        const double ay = seg.a.y - domain / 2.0;
        const double bx = seg.b.x - domain / 2.0;
        const double by = seg.b.y - domain / 2.0;
        maxRadiusError = std::max(maxRadiusError, std::fabs(std::sqrt(ax * ax + ay * ay) - radius));
        maxRadiusError = std::max(maxRadiusError, std::fabs(std::sqrt(bx * bx + by * by) - radius));
        const double segDx = bx - ax;
        const double segDy = by - ay;
        totalLength += std::sqrt(segDx * segDx + segDy * segDy);
    }

    const double kPi = 3.14159265358979323846;
    // Measured directly (not guessed): at this resolution, max radius
    // error is ~1.2e-4 (grid spacing dx is ~0.031, two orders of
    // magnitude larger) and total contour length is within ~8e-4 of
    // 2*pi*r -- both ordinary linear-interpolation error, comfortable
    // margins below.
    AETHER_CHECK(maxRadiusError < 5e-4);
    AETHER_CHECK(nearlyEqual(totalLength, 2.0 * kPi * radius, 5e-3));
}

// Marching cubes validated the same way as marching squares, one dimension
// up: f(x,y,z) = distance from a center point has the sphere of radius r
// as its f = r iso-surface, for any r. Every extracted triangle vertex
// must lie within ordinary grid-resolution error of that sphere, total
// surface area must approximate 4*pi*r^2, and -- the property specific to
// marching *tetrahedra*'s per-case orientation logic (see the class's own
// header comment) -- every triangle's normal must point away from the
// sphere's center, since the field increases outward.
void testMarchingCubesExtractsSphere() {
    const std::size_t n = 30;
    const double domain = 2.0; // covers [0, 2] in x, y, and z
    const double dx = domain / static_cast<double>(n);
    const double radius = 0.6;
    const Vector3 center(1.0, 1.0, 1.0);

    std::vector<double> field(n * n * n);
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double x = (static_cast<double>(i) + 0.5) * dx;
                const double y = (static_cast<double>(j) + 0.5) * dx;
                const double z = (static_cast<double>(k) + 0.5) * dx;
                const double dxc = x - center.x;
                const double dyc = y - center.y;
                const double dzc = z - center.z;
                field[i + j * n + k * n * n] = std::sqrt(dxc * dxc + dyc * dyc + dzc * dzc);
            }
        }
    }

    const auto triangles = marchingCubes3D(n, n, n, domain, domain, domain, field, radius);
    AETHER_CHECK(!triangles.empty());

    double maxRadiusError = 0.0;
    double totalArea = 0.0;
    for (const Triangle3D& tri : triangles) {
        for (const Vector3& p : {tri.a, tri.b, tri.c}) {
            const double d = (p - center).norm();
            maxRadiusError = std::max(maxRadiusError, std::fabs(d - radius));
        }
        const Vector3 normal = (tri.b - tri.a).cross(tri.c - tri.a);
        totalArea += 0.5 * normal.norm();

        const Vector3 centroid = (tri.a + tri.b + tri.c) * (1.0 / 3.0);
        // Outward orientation: the field increases away from `center`, so
        // every triangle's normal must point away from it too.
        AETHER_CHECK(normal.dot(centroid - center) > 0.0);
    }

    const double kPi = 3.14159265358979323846;
    // Measured directly (not guessed): at this resolution, max radius
    // error is ~2.6e-3 (grid spacing dx is ~0.067, more than an order of
    // magnitude larger) and total surface area is within ~0.02 of
    // 4*pi*r^2 (~4.524) -- both ordinary linear-interpolation/faceting
    // error, comfortable margins below.
    AETHER_CHECK(maxRadiusError < 1e-2);
    AETHER_CHECK(nearlyEqual(totalArea, 4.0 * kPi * radius * radius, 0.05));
}

// **A file format is only exported correctly if it can be read back and
// still say the same thing**, so this parses its own output rather than
// checking that the writer ran without throwing. There is no external VTK
// library here to validate against, and pulling one in to test a few
// hundred lines of text output would be a much larger dependency decision
// than the feature warrants -- but the round-trip is self-contained and
// decisive on its own: every count in the header, every coordinate, every
// connectivity index and every field value has to match the mesh the
// writer was given.
void testVtkWriterRoundTripsMeshAndFields() {
    // Two tetrahedra sharing a face, built by hand rather than generated:
    // the assertions below compare against literal expected values, which
    // only means something if the input is known exactly.
    DelaunayTetrahedralization3D tets;
    tets.addPoint(0.0, 0.0, 0.0);
    tets.addPoint(1.0, 0.0, 0.0);
    tets.addPoint(0.0, 1.0, 0.0);
    tets.addPoint(0.0, 0.0, 1.0);
    tets.addPoint(0.75, 0.75, 0.75);
    tets.tetrahedralize();
    const TetrahedralMesh mesh = TetrahedralMesh::fromTetrahedralization(tets);
    AETHER_CHECK(mesh.cellCount() > 0);

    // Values chosen so a transposed or off-by-one write shows up: each
    // cell's scalar is its own index, and each component of its vector is
    // a different function of that index, so no permutation of them
    // reproduces the original.
    CellScalarField pressure{"pressure", {}};
    CellVectorField velocity{"velocity", {}};
    for (std::size_t c = 0; c < mesh.cellCount(); ++c) {
        pressure.values.push_back(static_cast<double>(c) * 0.125);
        velocity.values.push_back(Vector3(static_cast<double>(c), -static_cast<double>(c) * 2.0,
                                           static_cast<double>(c) * 0.5 + 1.0));
    }

    const std::string path = "aether_vtk_roundtrip.vtk";
    writeTetrahedralMeshVtk(path, mesh, {pressure}, {velocity});

    std::ifstream in(path);
    AETHER_CHECK(in.good());
    std::string token;

    // Header: skip the four fixed lines, then read the POINTS block.
    std::string line;
    for (int i = 0; i < 4; ++i) {
        std::getline(in, line);
    }
    std::size_t pointCount = 0;
    in >> token >> pointCount >> token; // "POINTS" <n> "double"
    AETHER_CHECK(token == "double");
    AETHER_CHECK(pointCount == mesh.vertexCount());
    for (std::size_t v = 0; v < pointCount; ++v) {
        double x = 0.0, y = 0.0, z = 0.0;
        in >> x >> y >> z;
        const Vector3& expected = mesh.vertex(v);
        // Exact equality, not a tolerance: 17 significant digits round-trip
        // a double bit for bit, so anything else is a writer bug.
        AETHER_CHECK(x == expected.x && y == expected.y && z == expected.z);
    }

    std::size_t cellCount = 0;
    std::size_t totalInts = 0;
    in >> token >> cellCount >> totalInts; // "CELLS" <n> <n*5>
    AETHER_CHECK(cellCount == mesh.cellCount());
    AETHER_CHECK(totalInts == cellCount * 5);
    for (std::size_t c = 0; c < cellCount; ++c) {
        std::size_t perCell = 0;
        in >> perCell;
        AETHER_CHECK(perCell == 4);
        const std::vector<std::size_t>& expected = mesh.cellVertices(c);
        for (std::size_t i = 0; i < 4; ++i) {
            std::size_t index = 0;
            in >> index;
            AETHER_CHECK(index == expected[i]);
            AETHER_CHECK(index < pointCount);
        }
    }

    std::size_t typeCount = 0;
    in >> token >> typeCount; // "CELL_TYPES" <n>
    AETHER_CHECK(typeCount == cellCount);
    for (std::size_t c = 0; c < typeCount; ++c) {
        int type = 0;
        in >> type;
        AETHER_CHECK(type == 10); // VTK_TETRA
    }

    std::size_t dataCount = 0;
    in >> token >> dataCount; // "CELL_DATA" <n>
    AETHER_CHECK(dataCount == cellCount);

    std::string name, kind;
    int components = 0;
    in >> token >> name >> kind >> components; // "SCALARS" <name> "double" 1
    AETHER_CHECK(token == "SCALARS" && name == "pressure" && components == 1);
    std::getline(in, line);
    std::getline(in, line); // "LOOKUP_TABLE default"
    for (std::size_t c = 0; c < cellCount; ++c) {
        double value = 0.0;
        in >> value;
        AETHER_CHECK(value == pressure.values[c]);
    }

    in >> token >> name >> kind; // "VECTORS" <name> "double"
    AETHER_CHECK(token == "VECTORS" && name == "velocity");
    for (std::size_t c = 0; c < cellCount; ++c) {
        double x = 0.0, y = 0.0, z = 0.0;
        in >> x >> y >> z;
        AETHER_CHECK(x == velocity.values[c].x && y == velocity.values[c].y &&
                      z == velocity.values[c].z);
    }

    in.close();
    std::remove(path.c_str());

    std::printf("  [postprocessing_tests] VTK: %zu pontos e %zu celulas com 2 campos "
                "conferem apos releitura\n",
                pointCount, cellCount);
    std::fflush(stdout);
}

// A mismatched field length or a name with a space in it produces a file
// that another tool misparses rather than rejects, so both are refused at
// the call instead. Checked because a validation nobody exercises is a
// validation that quietly stops working.
void testVtkWriterRefusesInconsistentInput() {
    DelaunayTetrahedralization3D tets;
    tets.addPoint(0.0, 0.0, 0.0);
    tets.addPoint(1.0, 0.0, 0.0);
    tets.addPoint(0.0, 1.0, 0.0);
    tets.addPoint(0.0, 0.0, 1.0);
    tets.tetrahedralize();
    const TetrahedralMesh mesh = TetrahedralMesh::fromTetrahedralization(tets);

    bool refusedShortField = false;
    try {
        writeTetrahedralMeshVtk("aether_vtk_should_not_exist.vtk", mesh,
                                 {CellScalarField{"p", std::vector<double>(mesh.cellCount() + 1, 0.0)}});
    } catch (const std::invalid_argument&) {
        refusedShortField = true;
    }
    AETHER_CHECK(refusedShortField);

    bool refusedBadName = false;
    try {
        writeTetrahedralMeshVtk("aether_vtk_should_not_exist.vtk", mesh,
                                 {CellScalarField{"bad name",
                                                   std::vector<double>(mesh.cellCount(), 0.0)}});
    } catch (const std::invalid_argument&) {
        refusedBadName = true;
    }
    AETHER_CHECK(refusedBadName);

    std::printf("  [postprocessing_tests] VTK: campo de tamanho errado e nome com espaco recusados\n");
    std::fflush(stdout);
}

} // namespace

int main() {
    testStreamlineConservesTaylorGreenStreamFunction();
    testMarchingSquaresExtractsCircle();
    testMarchingCubesExtractsSphere();
    testVtkWriterRoundTripsMeshAndFields();
    testVtkWriterRefusesInconsistentInput();
    std::puts("aether_postprocessing_tests: all tests passed");
    return 0;
}
