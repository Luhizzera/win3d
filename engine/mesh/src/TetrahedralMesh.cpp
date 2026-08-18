#include "aether/mesh/TetrahedralMesh.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace aether::mesh {

using core::Vector3;

namespace {

// The four triangular faces of a tetrahedron (v0,v1,v2,v3), each wound so
// that (p1-p0) x (p2-p0) points *out* of the tetrahedron when it is
// positively oriented, i.e. when
// det = (v1-v0) . ((v2-v0) x (v3-v0)) > 0.
//
// Derivation for the last row, since an inverted table is the easiest
// mistake to make here and the hardest to notice: for the face opposite
// v3, wound (v0,v2,v1), the normal is
//   n = (v2-v0) x (v1-v0) = -[(v1-v0) x (v2-v0)]
// and outward means n points away from the remaining vertex v3:
//   n . (v3-v0) = -[(v1-v0) x (v2-v0)] . (v3-v0) = -det < 0   (outward)
// The naive winding (v0,v1,v2) gives exactly +det > 0, i.e. INWARD. The
// other three rows follow the same way.
//
// **This table was wrong in its first version -- every row reversed, so
// every normal pointed inward.** Neither the closure identity nor the
// owner/neighbour antisymmetry check catches that: a globally flipped set
// of normals still sums to zero over a cell, and still negates correctly
// between the two sides of a face. Only
// testTetrahedralMeshAreaVectorsPointOwnerToNeighbour() catches it, which
// is precisely why that test exists as a separate check.
constexpr std::array<std::array<int, 3>, 4> kFaceVertexTable = {{
    {{1, 2, 3}}, // opposite vertex 0
    {{0, 3, 2}}, // opposite vertex 1
    {{0, 1, 3}}, // opposite vertex 2
    {{0, 2, 1}}, // opposite vertex 3
}};

std::array<std::size_t, 3> sortedKey(const std::array<std::size_t, 3>& v) {
    std::array<std::size_t, 3> key = v;
    std::sort(key.begin(), key.end());
    return key;
}

} // namespace

TetrahedralMesh TetrahedralMesh::fromTetrahedralization(
    const DelaunayTetrahedralization3D& tetrahedralization) {
    TetrahedralMesh mesh;

    for (std::size_t i = 0; i < tetrahedralization.pointCount(); ++i) {
        mesh.mesh_.addVertex(tetrahedralization.point(i));
    }

    // First pass: cell geometry. A tetrahedron's centroid is exactly the
    // average of its vertices and its volume exactly |det|/6, so neither
    // is an approximation here.
    std::vector<std::array<std::size_t, 4>> cellVertices;
    for (std::size_t t = 0; t < tetrahedralization.tetrahedronCount(); ++t) {
        const auto& tet = tetrahedralization.tetrahedron(t);
        const Vector3& a = mesh.mesh_.vertex(tet.vertices[0]);
        const Vector3& b = mesh.mesh_.vertex(tet.vertices[1]);
        const Vector3& c = mesh.mesh_.vertex(tet.vertices[2]);
        const Vector3& d = mesh.mesh_.vertex(tet.vertices[3]);

        const double signedVolume = (b - a).dot((c - a).cross(d - a)) / 6.0;
        if (std::fabs(signedVolume) == 0.0) {
            continue; // degenerate: carries no flux, and would divide by zero downstream
        }

        cellVertices.push_back(tet.vertices);
        // Registered in the core mesh too, so a ScalarField built over
        // coreMesh() is indexed by exactly these cells, in this order.
        mesh.mesh_.addCell({tet.vertices[0], tet.vertices[1], tet.vertices[2], tet.vertices[3]});
        mesh.cellVolumes_.push_back(std::fabs(signedVolume));
        mesh.cellCentroids_.push_back((a + b + c + d) / 4.0);
    }
    mesh.cellFaces_.resize(mesh.cellVolumes_.size());

    // Second pass: faces. A face is identified by its sorted vertex triple,
    // so the two tetrahedra sharing it -- which list it in opposite winding
    // orders -- map to the same key.
    std::map<std::array<std::size_t, 3>, std::size_t> faceIndexByKey;

    for (std::size_t cell = 0; cell < cellVertices.size(); ++cell) {
        const auto& tv = cellVertices[cell];

        // A tetrahedron from this generator is positively oriented by
        // construction, but a mesh assembled some other way might not be;
        // checking costs one determinant and makes the winding table below
        // correct either way.
        const Vector3& a = mesh.mesh_.vertex(tv[0]);
        const Vector3& b = mesh.mesh_.vertex(tv[1]);
        const Vector3& c = mesh.mesh_.vertex(tv[2]);
        const Vector3& d = mesh.mesh_.vertex(tv[3]);
        const bool positivelyOriented = (b - a).dot((c - a).cross(d - a)) > 0.0;

        for (const auto& localFace : kFaceVertexTable) {
            std::array<std::size_t, 3> faceVerts = {tv[localFace[0]], tv[localFace[1]], tv[localFace[2]]};
            if (!positivelyOriented) {
                std::swap(faceVerts[1], faceVerts[2]);
            }

            const auto key = sortedKey(faceVerts);
            const auto existing = faceIndexByKey.find(key);
            if (existing != faceIndexByKey.end()) {
                // Second (and final) tetrahedron on this face: it becomes
                // the neighbour. The area vector already points away from
                // the owner, which is now exactly "towards the neighbour".
                Face& face = mesh.faces_[existing->second];
                face.neighbour = cell;
                mesh.cellFaces_[cell].push_back(existing->second);
                continue;
            }

            const Vector3& p0 = mesh.mesh_.vertex(faceVerts[0]);
            const Vector3& p1 = mesh.mesh_.vertex(faceVerts[1]);
            const Vector3& p2 = mesh.mesh_.vertex(faceVerts[2]);

            Face face;
            face.vertices = faceVerts;
            face.owner = cell;
            face.neighbour = kNoNeighbour;
            // Half the cross product is the triangle's area times its unit
            // normal, in one step -- no normalize-then-rescale, so the
            // divergence-theorem identity stays exact to roundoff.
            face.areaVector = (p1 - p0).cross(p2 - p0) * 0.5;
            face.centroid = (p0 + p1 + p2) / 3.0;

            faceIndexByKey.emplace(key, mesh.faces_.size());
            mesh.cellFaces_[cell].push_back(mesh.faces_.size());
            mesh.faces_.push_back(face);
        }
    }

    return mesh;
}

std::size_t TetrahedralMesh::boundaryFaceCount() const {
    std::size_t count = 0;
    for (const Face& face : faces_) {
        if (face.neighbour == kNoNeighbour) {
            ++count;
        }
    }
    return count;
}

Vector3 TetrahedralMesh::outwardAreaVector(std::size_t cell, std::size_t faceIndex) const {
    const Face& face = faces_.at(faceIndex);
    return face.owner == cell ? face.areaVector : -face.areaVector;
}

Vector3 TetrahedralMesh::cellAreaVectorSum(std::size_t cell) const {
    Vector3 sum;
    for (std::size_t faceIndex : cellFaces_.at(cell)) {
        sum += outwardAreaVector(cell, faceIndex);
    }
    return sum;
}

double TetrahedralMesh::totalVolume() const {
    double total = 0.0;
    for (double volume : cellVolumes_) {
        total += volume;
    }
    return total;
}

} // namespace aether::mesh
