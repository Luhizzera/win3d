#include "aether/postprocessing/MarchingCubes3D.hpp"
#include "aether/postprocessing/MarchingSquares2D.hpp"
#include "aether/postprocessing/Streamline2D.hpp"
#include "aether/testing/Check.hpp"

#include <cmath>
#include <cstdio>

using namespace aether::core;
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

} // namespace

int main() {
    testStreamlineConservesTaylorGreenStreamFunction();
    testMarchingSquaresExtractsCircle();
    testMarchingCubesExtractsSphere();
    std::puts("aether_postprocessing_tests: all tests passed");
    return 0;
}
