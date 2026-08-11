#include "aether/geometry/TriangleMesh.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <unordered_map>

namespace aether::geometry {

using aether::core::Vector3;

std::size_t TriangleMesh::addVertex(const Vector3& position) {
    vertices_.push_back(position);
    return vertices_.size() - 1;
}

std::size_t TriangleMesh::addTriangle(std::size_t a, std::size_t b, std::size_t c) {
    if (a >= vertices_.size() || b >= vertices_.size() || c >= vertices_.size()) {
        throw std::out_of_range("TriangleMesh::addTriangle: vertex index out of range");
    }
    triangles_.push_back(Triangle{{a, b, c}});
    return triangles_.size() - 1;
}

Vector3 TriangleMesh::faceNormal(std::size_t triangleIndex) const {
    const Triangle& t = triangle(triangleIndex);
    const Vector3& v0 = vertices_[t.vertices[0]];
    const Vector3& v1 = vertices_[t.vertices[1]];
    const Vector3& v2 = vertices_[t.vertices[2]];
    return (v1 - v0).cross(v2 - v0).normalized();
}

double TriangleMesh::faceArea(std::size_t triangleIndex) const {
    const Triangle& t = triangle(triangleIndex);
    const Vector3& v0 = vertices_[t.vertices[0]];
    const Vector3& v1 = vertices_[t.vertices[1]];
    const Vector3& v2 = vertices_[t.vertices[2]];
    return 0.5 * (v1 - v0).cross(v2 - v0).norm();
}

double TriangleMesh::surfaceArea() const {
    double total = 0.0;
    for (std::size_t i = 0; i < triangles_.size(); ++i) {
        total += faceArea(i);
    }
    return total;
}

double TriangleMesh::volume() const {
    double sixVolume = 0.0;
    for (const Triangle& t : triangles_) {
        const Vector3& v0 = vertices_[t.vertices[0]];
        const Vector3& v1 = vertices_[t.vertices[1]];
        const Vector3& v2 = vertices_[t.vertices[2]];
        sixVolume += v0.dot(v1.cross(v2));
    }
    return std::fabs(sixVolume) / 6.0;
}

std::size_t TriangleMesh::weldVertices(double tolerance) {
    struct Key {
        long long x, y, z;
        bool operator==(const Key& o) const { return x == o.x && y == o.y && z == o.z; }
    };
    struct KeyHash {
        std::size_t operator()(const Key& k) const {
            std::size_t h = std::hash<long long>{}(k.x);
            h ^= std::hash<long long>{}(k.y) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            h ^= std::hash<long long>{}(k.z) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
            return h;
        }
    };

    const double invTolerance = 1.0 / tolerance;
    auto quantize = [invTolerance](double v) { return std::llround(v * invTolerance); };

    std::unordered_map<Key, std::size_t, KeyHash> canonical;
    std::vector<std::size_t> remap(vertices_.size());
    std::vector<Vector3> newVertices;
    newVertices.reserve(vertices_.size());

    for (std::size_t i = 0; i < vertices_.size(); ++i) {
        const Vector3& v = vertices_[i];
        Key key{quantize(v.x), quantize(v.y), quantize(v.z)};
        auto it = canonical.find(key);
        if (it == canonical.end()) {
            const std::size_t newIndex = newVertices.size();
            newVertices.push_back(v);
            canonical.emplace(key, newIndex);
            remap[i] = newIndex;
        } else {
            remap[i] = it->second;
        }
    }

    const std::size_t merged = vertices_.size() - newVertices.size();
    vertices_ = std::move(newVertices);
    for (Triangle& t : triangles_) {
        t.vertices[0] = remap[t.vertices[0]];
        t.vertices[1] = remap[t.vertices[1]];
        t.vertices[2] = remap[t.vertices[2]];
    }
    return merged;
}

std::size_t TriangleMesh::removeDegenerateTriangles(double areaTolerance) {
    const std::size_t before = triangles_.size();
    std::vector<Triangle> kept;
    kept.reserve(triangles_.size());
    for (std::size_t i = 0; i < triangles_.size(); ++i) {
        const Triangle& t = triangles_[i];
        const bool repeatedVertex =
            t.vertices[0] == t.vertices[1] || t.vertices[1] == t.vertices[2] || t.vertices[0] == t.vertices[2];
        if (!repeatedVertex && faceArea(i) > areaTolerance) {
            kept.push_back(t);
        }
    }
    triangles_ = std::move(kept);
    return before - triangles_.size();
}

void TriangleMesh::reorientNormalsOutward() {
    if (vertices_.empty()) {
        return;
    }
    Vector3 centroid;
    for (const Vector3& v : vertices_) {
        centroid += v;
    }
    centroid = centroid / static_cast<double>(vertices_.size());

    for (Triangle& t : triangles_) {
        const Vector3& v0 = vertices_[t.vertices[0]];
        const Vector3& v1 = vertices_[t.vertices[1]];
        const Vector3& v2 = vertices_[t.vertices[2]];
        const Vector3 faceCentroid = (v0 + v1 + v2) / 3.0;
        const Vector3 normal = (v1 - v0).cross(v2 - v0);
        if (normal.dot(faceCentroid - centroid) < 0.0) {
            std::swap(t.vertices[1], t.vertices[2]);
        }
    }
}

std::vector<TriangleMesh::Edge> TriangleMesh::findBoundaryEdges() const {
    std::map<std::pair<std::size_t, std::size_t>, int> edgeCount;
    auto addEdge = [&edgeCount](std::size_t a, std::size_t b) {
        if (a > b) {
            std::swap(a, b);
        }
        edgeCount[{a, b}] += 1;
    };

    for (const Triangle& t : triangles_) {
        addEdge(t.vertices[0], t.vertices[1]);
        addEdge(t.vertices[1], t.vertices[2]);
        addEdge(t.vertices[2], t.vertices[0]);
    }

    std::vector<Edge> boundary;
    for (const auto& [edge, count] : edgeCount) {
        if (count == 1) {
            boundary.push_back(Edge{edge.first, edge.second});
        }
    }
    return boundary;
}

} // namespace aether::geometry
