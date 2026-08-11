#pragma once

#include "aether/core/Vector3.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace aether::mesh {

// Module 3, first genuinely unstructured mesh: 2D Delaunay triangulation
// of a point set via the Bowyer-Watson incremental algorithm. This fills a
// gap left open since Module 3's first pass (StructuredGrid3D only covers
// uniform Cartesian grids); a general 3D tetrahedralizer and boundary-
// constrained triangulation (respecting an imported geometry's edges, not
// just its point cloud) are natural next steps built on this same
// algorithmic core.
//
// Points lie in the z=0 plane (stored as Vector3 for compatibility with
// the rest of the engine, z always 0). Triangles are always wound
// counter-clockwise.
class DelaunayTriangulation2D {
public:
    struct Triangle {
        std::array<std::size_t, 3> vertices;
    };

    // Adds a point (before triangulate() is called) and returns its index.
    std::size_t addPoint(double x, double y);

    // Computes the Delaunay triangulation of every point added so far,
    // via Bowyer-Watson: insert points one at a time, remove every
    // triangle whose circumcircle contains the new point (the resulting
    // hole is always star-shaped around the new point for a valid
    // triangulation), and re-triangulate the hole by connecting the new
    // point to its boundary. A large enclosing "super triangle" seeds the
    // process and is discarded at the end.
    void triangulate();

    std::size_t pointCount() const { return points_.size(); }
    std::size_t triangleCount() const { return triangles_.size(); }
    const core::Vector3& point(std::size_t i) const { return points_.at(i); }
    const Triangle& triangle(std::size_t i) const { return triangles_.at(i); }

    // Brute-force check of the defining Delaunay property: no point lies
    // strictly inside any triangle's circumcircle. O(n * triangleCount);
    // intended for validating modest point sets in tests, not as a
    // runtime query.
    bool satisfiesDelaunayProperty(double tolerance = 1e-7) const;

private:
    double inCircumcircle(const Triangle& t, std::size_t pointIndex) const;

    std::vector<core::Vector3> points_;
    std::vector<Triangle> triangles_;
};

} // namespace aether::mesh
