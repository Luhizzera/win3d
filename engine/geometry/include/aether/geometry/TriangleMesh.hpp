#pragma once

#include "aether/core/Vector3.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace aether::geometry {

struct Triangle {
    std::array<std::size_t, 3> vertices;
};

// Indexed triangle surface mesh (e.g. loaded from STL) plus the small set of
// repair/inspection operations Module 2 of the plan calls for: normal
// correction, hole detection, degenerate-triangle removal, volume/area.
// This is the *geometry* layer; volumetric mesh generation (tet/hex/prism,
// AMR) is a separate, later module and does not belong here.
class TriangleMesh {
public:
    std::size_t addVertex(const core::Vector3& position);
    std::size_t addTriangle(std::size_t a, std::size_t b, std::size_t c);

    std::size_t vertexCount() const { return vertices_.size(); }
    std::size_t triangleCount() const { return triangles_.size(); }

    const core::Vector3& vertex(std::size_t index) const { return vertices_.at(index); }
    const Triangle& triangle(std::size_t index) const { return triangles_.at(index); }

    core::Vector3 faceNormal(std::size_t triangleIndex) const;
    double faceArea(std::size_t triangleIndex) const;

    double surfaceArea() const;

    // Unsigned volume via the divergence theorem (sum of signed tetrahedra
    // from the origin); does not require prior normal correction.
    double volume() const;

    // Merges vertices within `tolerance` of each other and remaps triangle
    // indices accordingly. STL files store 3 independent vertices per
    // triangle with no shared connectivity, so this must run before
    // findBoundaryEdges()/isWatertight() can detect anything. Returns the
    // number of vertices merged away.
    std::size_t weldVertices(double tolerance = 1e-6);

    // Removes triangles that are degenerate (repeated vertex index, or face
    // area below areaTolerance). Returns how many were removed.
    std::size_t removeDegenerateTriangles(double areaTolerance = 1e-12);

    // Flips triangle winding so each face normal points away from the
    // mesh's vertex centroid. A simple heuristic that is correct for
    // star-shaped/mostly-convex closed meshes; it is not a general solution
    // for arbitrary non-convex geometry.
    void reorientNormalsOutward();

    struct Edge {
        std::size_t a;
        std::size_t b;
    };

    // Edges referenced by exactly one triangle (requires weldVertices() to
    // have been run for shared edges to be recognized as such). A closed,
    // watertight mesh has none.
    std::vector<Edge> findBoundaryEdges() const;
    bool isWatertight() const { return findBoundaryEdges().empty(); }

private:
    std::vector<core::Vector3> vertices_;
    std::vector<Triangle> triangles_;
};

} // namespace aether::geometry
