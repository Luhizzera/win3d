#pragma once

#include "aether/core/Vector3.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace aether::mesh {

// Module 3: triangulates a simple (non-self-intersecting), CCW-wound
// polygon by ear clipping followed by Delaunay edge flips restricted to
// never touch a boundary edge (Lawson's algorithm, constrained) -- the
// standard combination that yields a genuine Constrained Delaunay
// Triangulation (CDT): every input boundary edge survives exactly (what
// DelaunayTriangulation2D deliberately does not provide, since it
// triangulates a raw point cloud and may not preserve any particular
// edge), *and* the interior diagonals are chosen for Delaunay quality
// (unlike plain ear clipping alone, which can leave thin/sliver
// triangles). This is the piece actually needed to turn an imported 2D
// boundary (e.g. a cross-section of TriangleMesh/STL geometry) into a
// mesh a solver can use with reasonable triangle quality.
//
// Ear clipping repeatedly finds a convex vertex whose triangle with its
// two polygon neighbors contains no other polygon vertex (an "ear"),
// emits that triangle, and removes the vertex from the polygon, until
// only a single triangle remains. O(n^2) worst case (each ear search scans
// the remaining vertices) -- fine for the modest boundaries this is meant
// for; not the O(n log n) algorithm real production meshers use. The
// subsequent edge-flip pass repeatedly finds a non-boundary edge shared by
// two triangles that violates the Delaunay in-circle criterion and flips
// it, until none remain (or a safety iteration cap is hit).
//
// Does not (yet) handle holes/multiple contours, and does not insert new
// (Steiner) points -- so triangle quality is still limited by whatever
// diagonals are possible using only the input vertices. A full CDT
// refinement (Steiner point insertion for guaranteed angle/size bounds)
// and a 3D tetrahedralizer are natural next steps on this same core.
class PolygonTriangulation2D {
public:
    struct Triangle {
        std::array<std::size_t, 3> vertices;
    };

    // Adds a boundary vertex in order; the polygon is implicitly closed
    // (the last vertex connects back to the first).
    std::size_t addVertex(double x, double y);

    // Triangulates via ear clipping, then improves the result to a
    // Constrained Delaunay Triangulation via boundary-preserving edge
    // flips. The polygon must be simple and wound counter-clockwise.
    void triangulate();

    std::size_t vertexCount() const { return vertices_.size(); }
    std::size_t triangleCount() const { return triangles_.size(); }
    const core::Vector3& vertex(std::size_t i) const { return vertices_.at(i); }
    const Triangle& triangle(std::size_t i) const { return triangles_.at(i); }

    // Signed area of the input polygon (shoelace formula); positive for a
    // CCW-wound simple polygon.
    double polygonArea() const;

    // True if every *non-boundary* (diagonal) edge satisfies the Delaunay
    // in-circle criterion -- boundary edges are intentionally exempt,
    // since a CDT is only required to be locally Delaunay away from its
    // constrained edges.
    bool isLocallyDelaunay(double tolerance = 1e-7) const;

private:
    void flipToLocalDelaunay();
    bool isBoundaryEdge(std::size_t a, std::size_t b) const;

    std::vector<core::Vector3> vertices_;
    std::vector<Triangle> triangles_;
};

} // namespace aether::mesh
