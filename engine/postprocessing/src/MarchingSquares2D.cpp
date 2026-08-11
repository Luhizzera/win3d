#include "aether/postprocessing/MarchingSquares2D.hpp"

namespace aether::postprocessing {

using aether::core::Vector3;

namespace {

Vector3 cellCenter(std::size_t i, std::size_t j, double dx, double dy) {
    return Vector3((static_cast<double>(i) + 0.5) * dx, (static_cast<double>(j) + 0.5) * dy, 0.0);
}

Vector3 interpolateCrossing(const Vector3& a, double valueA, const Vector3& b, double valueB,
                             double isoValue) {
    const double t = (isoValue - valueA) / (valueB - valueA);
    return a + (b - a) * t;
}

} // namespace

std::vector<Segment2D> marchingSquares2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                                          const std::vector<double>& field, double isoValue) {
    const double dx = lengthX / static_cast<double>(nx);
    const double dy = lengthY / static_cast<double>(ny);

    std::vector<Segment2D> segments;
    if (nx < 2 || ny < 2) {
        return segments;
    }

    for (std::size_t j = 0; j + 1 < ny; ++j) {
        for (std::size_t i = 0; i + 1 < nx; ++i) {
            const Vector3 c0 = cellCenter(i, j, dx, dy);
            const Vector3 c1 = cellCenter(i + 1, j, dx, dy);
            const Vector3 c2 = cellCenter(i + 1, j + 1, dx, dy);
            const Vector3 c3 = cellCenter(i, j + 1, dx, dy);
            const double v0 = field[i + j * nx];
            const double v1 = field[(i + 1) + j * nx];
            const double v2 = field[(i + 1) + (j + 1) * nx];
            const double v3 = field[i + (j + 1) * nx];

            const bool above[4] = {v0 > isoValue, v1 > isoValue, v2 > isoValue, v3 > isoValue};
            const bool edgeCrossed[4] = {above[0] != above[1], above[1] != above[2], above[2] != above[3],
                                          above[3] != above[0]};

            Vector3 edgePoint[4];
            int numCrossed = 0;
            if (edgeCrossed[0]) {
                edgePoint[0] = interpolateCrossing(c0, v0, c1, v1, isoValue);
                ++numCrossed;
            }
            if (edgeCrossed[1]) {
                edgePoint[1] = interpolateCrossing(c1, v1, c2, v2, isoValue);
                ++numCrossed;
            }
            if (edgeCrossed[2]) {
                edgePoint[2] = interpolateCrossing(c2, v2, c3, v3, isoValue);
                ++numCrossed;
            }
            if (edgeCrossed[3]) {
                edgePoint[3] = interpolateCrossing(c3, v3, c0, v0, isoValue);
                ++numCrossed;
            }

            if (numCrossed == 0) {
                continue;
            }
            if (numCrossed == 2) {
                // Exactly one pair of crossed edges -- connect them.
                int first = -1;
                for (int e = 0; e < 4; ++e) {
                    if (edgeCrossed[e]) {
                        if (first < 0) {
                            first = e;
                        } else {
                            segments.push_back({edgePoint[first], edgePoint[e]});
                        }
                    }
                }
                continue;
            }
            // numCrossed == 4: the ambiguous "saddle" case (diagonal
            // corners share state, e.g. corners 0,2 above and 1,3 below).
            // Resolved with the average-corner-value heuristic -- see the
            // header comment.
            const double average = (v0 + v1 + v2 + v3) / 4.0;
            if ((average > isoValue) == above[0]) {
                segments.push_back({edgePoint[3], edgePoint[0]}); // isolates corner 0
                segments.push_back({edgePoint[1], edgePoint[2]}); // isolates corner 2
            } else {
                segments.push_back({edgePoint[0], edgePoint[1]}); // isolates corner 1
                segments.push_back({edgePoint[2], edgePoint[3]}); // isolates corner 3
            }
        }
    }
    return segments;
}

} // namespace aether::postprocessing
