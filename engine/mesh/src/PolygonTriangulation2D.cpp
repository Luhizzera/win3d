#include "aether/mesh/PolygonTriangulation2D.hpp"

#include "aether/mesh/RobustPredicates.hpp"

#include <cmath>
#include <map>
#include <utility>

namespace aether::mesh {

using aether::core::Vector3;
using Triangle = PolygonTriangulation2D::Triangle;

namespace {

// Every geometric sign test in this file goes through
// RobustPredicates' exact arithmetic rather than plain double
// determinants -- see RobustPredicates.hpp. These are all *sign* tests
// (turn direction, point-in-triangle, in-circumcircle), so exactness is
// exactly what matters and no magnitude is ever needed.

// Turn direction at b, going a -> b -> c: positive means a left turn
// (convex, for a CCW-wound polygon). Note orientation2D(a,b,c) computes
// the cross product of (b-a) and (c-a) where this historically used
// (b-a) x (c-b) -- algebraically the same sign, since (c-b) = (c-a)-(b-a)
// and the (b-a) x (b-a) term vanishes.
int turnCross(const Vector3& a, const Vector3& b, const Vector3& c) { return orientation2D(a, b, c); }

// True if p lies inside or on the boundary of triangle (a,b,c), regardless
// of the triangle's winding.
bool pointInOrOnTriangle(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c) {
    const int d1 = orientation2D(p, a, b);
    const int d2 = orientation2D(p, b, c);
    const int d3 = orientation2D(p, c, a);
    const bool hasNeg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    const bool hasPos = (d1 > 0) || (d2 > 0) || (d3 > 0);
    return !(hasNeg && hasPos);
}

// Positive iff p lies strictly inside the circumcircle of the CCW-wound
// triangle (a,b,c). Invariant under cyclic rotation of (a,b,c), so it can
// be applied directly to a triangle's stored vertex order without first
// figuring out which vertex is "opposite" anything.
int inCircumcircleDet(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& p) {
    return inCircle2D(a, b, c, p);
}

// Makes a CCW-wound triangle from three point indices, given their
// coordinates.
Triangle makeCcwTriangle(std::size_t a, std::size_t b, std::size_t c, const std::vector<Vector3>& points) {
    if (orientation2D(points[a], points[b], points[c]) < 0) {
        return {{a, c, b}};
    }
    return {{a, b, c}};
}

// The vertex of triangle t that is not a or b (t must contain edge a-b).
std::size_t oppositeVertex(const Triangle& t, std::size_t a, std::size_t b) {
    for (std::size_t v : t.vertices) {
        if (v != a && v != b) {
            return v;
        }
    }
    return t.vertices[0]; // unreachable if t genuinely contains edge (a,b)
}

// Maps each undirected edge to the (1 or 2) triangles that contain it.
std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>>
buildEdgeTriangleMap(const std::vector<Triangle>& triangles) {
    std::map<std::pair<std::size_t, std::size_t>, std::vector<std::size_t>> edgeTriangles;
    for (std::size_t t = 0; t < triangles.size(); ++t) {
        const auto& v = triangles[t].vertices;
        for (int e = 0; e < 3; ++e) {
            std::size_t a = v[static_cast<std::size_t>(e)];
            std::size_t b = v[static_cast<std::size_t>((e + 1) % 3)];
            if (a > b) {
                std::swap(a, b);
            }
            edgeTriangles[{a, b}].push_back(t);
        }
    }
    return edgeTriangles;
}

} // namespace

std::size_t PolygonTriangulation2D::addVertex(double x, double y) {
    vertices_.emplace_back(x, y, 0.0);
    return vertices_.size() - 1;
}

double PolygonTriangulation2D::polygonArea() const {
    double sum = 0.0;
    const std::size_t n = vertices_.size();
    for (std::size_t i = 0; i < n; ++i) {
        const Vector3& p1 = vertices_[i];
        const Vector3& p2 = vertices_[(i + 1) % n];
        sum += p1.x * p2.y - p2.x * p1.y;
    }
    return sum / 2.0;
}

void PolygonTriangulation2D::triangulate() {
    triangles_.clear();
    const std::size_t n = vertices_.size();
    if (n < 3) {
        return;
    }

    std::vector<std::size_t> indices(n);
    for (std::size_t i = 0; i < n; ++i) {
        indices[i] = i;
    }

    // Safety cap: a valid simple polygon clips exactly n-3 ears before the
    // final triangle. If no ear is found within that many attempts (a
    // self-intersecting or degenerate input), stop rather than looping
    // forever.
    std::size_t guard = 0;
    const std::size_t guardLimit = n * n + 8;

    while (indices.size() > 3 && guard < guardLimit) {
        ++guard;
        const std::size_t m = indices.size();
        bool earFound = false;

        for (std::size_t i = 0; i < m; ++i) {
            const std::size_t prevIdx = indices[(i + m - 1) % m];
            const std::size_t currIdx = indices[i];
            const std::size_t nextIdx = indices[(i + 1) % m];

            const Vector3& a = vertices_[prevIdx];
            const Vector3& b = vertices_[currIdx];
            const Vector3& c = vertices_[nextIdx];

            if (turnCross(a, b, c) <= 0) {
                continue; // reflex (concave) vertex: cannot be an ear
            }

            bool containsOther = false;
            for (std::size_t idx : indices) {
                if (idx == prevIdx || idx == currIdx || idx == nextIdx) {
                    continue;
                }
                if (pointInOrOnTriangle(vertices_[idx], a, b, c)) {
                    containsOther = true;
                    break;
                }
            }
            if (containsOther) {
                continue;
            }

            triangles_.push_back({{prevIdx, currIdx, nextIdx}});
            indices.erase(indices.begin() + static_cast<std::ptrdiff_t>(i));
            earFound = true;
            break;
        }

        if (!earFound) {
            break; // degenerate/self-intersecting input; stop rather than spin
        }
    }

    if (indices.size() == 3) {
        triangles_.push_back({{indices[0], indices[1], indices[2]}});
    }

    flipToLocalDelaunay();
}

bool PolygonTriangulation2D::isBoundaryEdge(std::size_t a, std::size_t b) const {
    const std::size_t n = vertices_.size();
    if (a > b) {
        std::swap(a, b);
    }
    return (b == a + 1) || (a == 0 && b == n - 1);
}

void PolygonTriangulation2D::flipToLocalDelaunay() {
    std::size_t safety = 0;
    const std::size_t safetyLimit = triangles_.size() * triangles_.size() + 16;

    bool anyFlip = true;
    while (anyFlip && safety < safetyLimit) {
        ++safety;
        anyFlip = false;

        const auto edgeTriangles = buildEdgeTriangleMap(triangles_);
        for (const auto& [edge, tris] : edgeTriangles) {
            if (tris.size() != 2 || isBoundaryEdge(edge.first, edge.second)) {
                continue;
            }

            const std::size_t t1 = tris[0];
            const std::size_t t2 = tris[1];
            const std::size_t opp2 = oppositeVertex(triangles_[t2], edge.first, edge.second);
            const auto& v1 = triangles_[t1].vertices;

            if (inCircumcircleDet(vertices_[v1[0]], vertices_[v1[1]], vertices_[v1[2]], vertices_[opp2]) <=
                0) {
                continue; // already locally Delaunay
            }

            const std::size_t opp1 = oppositeVertex(triangles_[t1], edge.first, edge.second);
            triangles_[t1] = makeCcwTriangle(opp1, edge.first, opp2, vertices_);
            triangles_[t2] = makeCcwTriangle(opp1, opp2, edge.second, vertices_);
            anyFlip = true;
        }
    }
}

bool PolygonTriangulation2D::isLocallyDelaunay(double tolerance) const {
    const auto edgeTriangles = buildEdgeTriangleMap(triangles_);
    for (const auto& [edge, tris] : edgeTriangles) {
        if (tris.size() != 2 || isBoundaryEdge(edge.first, edge.second)) {
            continue;
        }
        const std::size_t t1 = tris[0];
        const std::size_t t2 = tris[1];
        const std::size_t opp2 = oppositeVertex(triangles_[t2], edge.first, edge.second);
        const auto& v1 = triangles_[t1].vertices;
        if (inCircumcircleDet(vertices_[v1[0]], vertices_[v1[1]], vertices_[v1[2]], vertices_[opp2]) >
            tolerance) {
            return false;
        }
    }
    return true;
}

} // namespace aether::mesh
