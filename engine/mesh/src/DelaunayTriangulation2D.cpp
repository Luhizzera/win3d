#include "aether/mesh/DelaunayTriangulation2D.hpp"

#include "aether/mesh/RobustPredicates.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace aether::mesh {

using aether::core::Vector3;

namespace {

// Makes a CCW-wound triangle from three point indices, given the actual
// point coordinates (needed to determine winding). Winding is decided by
// orientation2D's exact sign rather than a plain double cross product --
// see RobustPredicates.hpp.
DelaunayTriangulation2D::Triangle makeTriangle(std::size_t a, std::size_t b, std::size_t c,
                                                const std::vector<Vector3>& points) {
    if (orientation2D(points[a], points[b], points[c]) < 0) {
        return {{a, c, b}};
    }
    return {{a, b, c}};
}

bool triangleHasVertex(const DelaunayTriangulation2D::Triangle& t, std::size_t v) {
    return t.vertices[0] == v || t.vertices[1] == v || t.vertices[2] == v;
}

struct Edge {
    std::size_t a;
    std::size_t b;
};

bool sameEdge(const Edge& lhs, const Edge& rhs) {
    return (lhs.a == rhs.a && lhs.b == rhs.b) || (lhs.a == rhs.b && lhs.b == rhs.a);
}

} // namespace

std::size_t DelaunayTriangulation2D::addPoint(double x, double y) {
    points_.emplace_back(x, y, 0.0);
    return points_.size() - 1;
}

// Routed through inCircle2D's exact arithmetic rather than a plain double
// 3x3 determinant (see RobustPredicates.hpp). Still returns a double
// (exactly -1.0, 0.0 or 1.0) so every existing `> 0.0` / `> tolerance`
// call site keeps working unchanged -- tolerances below 1 only ever
// distinguish "sign is exactly +1" from "0 or -1", which is the intended
// check.
double DelaunayTriangulation2D::inCircumcircle(const Triangle& t, std::size_t pointIndex) const {
    return static_cast<double>(inCircle2D(points_[t.vertices[0]], points_[t.vertices[1]],
                                           points_[t.vertices[2]], points_[pointIndex]));
}

void DelaunayTriangulation2D::triangulate() {
    const std::size_t originalPointCount = points_.size();
    triangles_.clear();
    if (originalPointCount == 0) {
        return;
    }

    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max();
    double maxY = std::numeric_limits<double>::lowest();
    for (std::size_t i = 0; i < originalPointCount; ++i) {
        minX = std::min(minX, points_[i].x);
        maxX = std::max(maxX, points_[i].x);
        minY = std::min(minY, points_[i].y);
        maxY = std::max(maxY, points_[i].y);
    }
    const double dx = maxX - minX;
    const double dy = maxY - minY;
    const double deltaMax = std::max(dx, dy) + 1.0; // +1 guards against a degenerate single-point case
    const double midX = (minX + maxX) / 2.0;
    const double midY = (minY + maxY) / 2.0;

    // Super triangle; its vertices are appended after the real points so
    // their indices are known, and every triangle touching them is
    // discarded at the end.
    //
    // As in DelaunayTetrahedralization3D (see that class's much longer
    // note on the same point), the scale here must be far larger than
    // merely "big enough to enclose every input point". The super
    // vertices are a *finite* stand-in for points at infinity, so a
    // triangle built on real points plus a super vertex gets a finite
    // circumcircle where the true one is a half-plane; if they are only
    // moderately far out, that circumcircle can wrongly test as
    // containing a nearby real point and delete a triangle that should
    // have survived. In 3D this produced a real, reproducible bug (a
    // missing hull-adjacent tetrahedron); the identical 20*deltaMax scale
    // used to be here too, and an A/B measurement confirmed the bug was
    // live here as well, not merely theoretical -- over 200 randomized
    // point sets, comparing the triangulated area against the
    // independently-computed convex hull area (a triangulation must
    // exactly fill its hull), the old scale failed 10 times (190/200)
    // while the scale below passes all 200. Cost is negligible: three
    // points, discarded at the end.
    const double r = 1.0e8 * deltaMax;
    const std::size_t superA = points_.size();
    points_.emplace_back(midX - r, midY - deltaMax, 0.0);
    const std::size_t superB = points_.size();
    points_.emplace_back(midX, midY + r, 0.0);
    const std::size_t superC = points_.size();
    points_.emplace_back(midX + r, midY - deltaMax, 0.0);

    triangles_.push_back(makeTriangle(superA, superB, superC, points_));

    for (std::size_t p = 0; p < originalPointCount; ++p) {
        std::vector<std::size_t> badTriangleIndices;
        for (std::size_t t = 0; t < triangles_.size(); ++t) {
            if (inCircumcircle(triangles_[t], p) > 0.0) {
                badTriangleIndices.push_back(t);
            }
        }

        // The hole left by removing the bad triangles is star-shaped
        // around p; its boundary is exactly the edges that belong to only
        // one bad triangle (edges shared between two bad triangles are
        // interior to the hole and disappear).
        std::vector<Edge> boundary;
        for (std::size_t bi : badTriangleIndices) {
            const Triangle& tri = triangles_[bi];
            const Edge edges[3] = {{tri.vertices[0], tri.vertices[1]},
                                    {tri.vertices[1], tri.vertices[2]},
                                    {tri.vertices[2], tri.vertices[0]}};
            for (const Edge& e : edges) {
                std::size_t sharedCount = 0;
                for (std::size_t bj : badTriangleIndices) {
                    if (bj == bi) {
                        continue;
                    }
                    const Triangle& other = triangles_[bj];
                    const Edge otherEdges[3] = {{other.vertices[0], other.vertices[1]},
                                                 {other.vertices[1], other.vertices[2]},
                                                 {other.vertices[2], other.vertices[0]}};
                    for (const Edge& oe : otherEdges) {
                        if (sameEdge(e, oe)) {
                            ++sharedCount;
                        }
                    }
                }
                if (sharedCount == 0) {
                    boundary.push_back(e);
                }
            }
        }

        std::vector<Triangle> kept;
        kept.reserve(triangles_.size());
        for (std::size_t t = 0; t < triangles_.size(); ++t) {
            if (std::find(badTriangleIndices.begin(), badTriangleIndices.end(), t) ==
                badTriangleIndices.end()) {
                kept.push_back(triangles_[t]);
            }
        }
        for (const Edge& e : boundary) {
            kept.push_back(makeTriangle(e.a, e.b, p, points_));
        }
        triangles_ = std::move(kept);
    }

    std::vector<Triangle> withoutSuper;
    withoutSuper.reserve(triangles_.size());
    for (const Triangle& t : triangles_) {
        if (!triangleHasVertex(t, superA) && !triangleHasVertex(t, superB) &&
            !triangleHasVertex(t, superC)) {
            withoutSuper.push_back(t);
        }
    }
    triangles_ = std::move(withoutSuper);

    points_.resize(originalPointCount); // drop the super triangle's vertices
}

bool DelaunayTriangulation2D::satisfiesDelaunayProperty(double tolerance) const {
    for (const Triangle& t : triangles_) {
        for (std::size_t p = 0; p < points_.size(); ++p) {
            if (triangleHasVertex(t, p)) {
                continue;
            }
            if (inCircumcircle(t, p) > tolerance) {
                return false;
            }
        }
    }
    return true;
}

} // namespace aether::mesh
