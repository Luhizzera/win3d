#include "aether/postprocessing/MarchingCubes3D.hpp"

#include <array>
#include <cstddef>

namespace aether::postprocessing {

using aether::core::Vector3;

namespace {

Vector3 cellCenter(std::size_t i, std::size_t j, std::size_t k, double dx, double dy, double dz) {
    return Vector3((static_cast<double>(i) + 0.5) * dx, (static_cast<double>(j) + 0.5) * dy,
                    (static_cast<double>(k) + 0.5) * dz);
}

Vector3 interpolateCrossing(const Vector3& a, double valueA, const Vector3& b, double valueB,
                             double isoValue) {
    const double t = (isoValue - valueA) / (valueB - valueA);
    return a + (b - a) * t;
}

// Appends the triangle (p0,p1,p2), flipping its winding if needed so its
// normal points toward `outwardTarget` (see the header comment: normal
// points from below isoValue toward above isoValue).
void pushOrientedTriangle(std::vector<Triangle3D>& triangles, const Vector3& p0, const Vector3& p1,
                           const Vector3& p2, const Vector3& outwardTarget) {
    const Vector3 centroid = (p0 + p1 + p2) * (1.0 / 3.0);
    const Vector3 normal = (p1 - p0).cross(p2 - p0);
    const double alignment = normal.dot(outwardTarget - centroid);
    if (alignment >= 0.0) {
        triangles.push_back({p0, p1, p2});
    } else {
        triangles.push_back({p0, p2, p1});
    }
}

// Processes a single tetrahedron with corners p[0..3] and field values
// v[0..3], appending 0, 1, or 2 triangles depending on how many corners
// are above isoValue -- see the header comment for the derivation.
void processTetrahedron(const std::array<Vector3, 4>& p, const std::array<double, 4>& v, double isoValue,
                         std::vector<Triangle3D>& triangles) {
    std::array<bool, 4> above{};
    int aboveCount = 0;
    for (int i = 0; i < 4; ++i) {
        above[i] = v[i] > isoValue;
        aboveCount += above[i] ? 1 : 0;
    }

    if (aboveCount == 0 || aboveCount == 4) {
        return;
    }

    if (aboveCount == 1 || aboveCount == 3) {
        int lone = -1;
        for (int i = 0; i < 4; ++i) {
            int matches = 0;
            for (int j = 0; j < 4; ++j) {
                if (above[j] == above[i]) {
                    ++matches;
                }
            }
            if (matches == 1) {
                lone = i;
                break;
            }
        }

        std::array<Vector3, 3> crossings;
        std::array<Vector3, 3> otherPoints;
        int count = 0;
        for (int i = 0; i < 4; ++i) {
            if (i == lone) {
                continue;
            }
            crossings[count] = interpolateCrossing(p[lone], v[lone], p[i], v[i], isoValue);
            otherPoints[count] = p[i];
            ++count;
        }

        const Vector3 outwardTarget =
            above[lone] ? p[lone] : (otherPoints[0] + otherPoints[1] + otherPoints[2]) * (1.0 / 3.0);
        pushOrientedTriangle(triangles, crossings[0], crossings[1], crossings[2], outwardTarget);
        return;
    }

    // aboveCount == 2: the "above" pair {a0,a1} and "below" pair {b0,b1}
    // form a quadrilateral via the 4 edges crossing between them, cyclic
    // order pA0B0 -> pA1B0 -> pA1B1 -> pA0B1 (see the header comment).
    std::array<int, 2> aboveIdx{};
    std::array<int, 2> belowIdx{};
    int aboveCursor = 0;
    int belowCursor = 0;
    for (int i = 0; i < 4; ++i) {
        if (above[i]) {
            aboveIdx[aboveCursor++] = i;
        } else {
            belowIdx[belowCursor++] = i;
        }
    }
    const int a0 = aboveIdx[0];
    const int a1 = aboveIdx[1];
    const int b0 = belowIdx[0];
    const int b1 = belowIdx[1];

    const Vector3 pA0B0 = interpolateCrossing(p[a0], v[a0], p[b0], v[b0], isoValue);
    const Vector3 pA1B0 = interpolateCrossing(p[a1], v[a1], p[b0], v[b0], isoValue);
    const Vector3 pA1B1 = interpolateCrossing(p[a1], v[a1], p[b1], v[b1], isoValue);
    const Vector3 pA0B1 = interpolateCrossing(p[a0], v[a0], p[b1], v[b1], isoValue);

    const Vector3 outwardTarget = (p[a0] + p[a1]) * 0.5;
    pushOrientedTriangle(triangles, pA0B0, pA1B0, pA1B1, outwardTarget);
    pushOrientedTriangle(triangles, pA0B0, pA1B1, pA0B1, outwardTarget);
}

} // namespace

std::vector<Triangle3D> marchingCubes3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                         double lengthY, double lengthZ, const std::vector<double>& field,
                                         double isoValue) {
    const double dx = lengthX / static_cast<double>(nx);
    const double dy = lengthY / static_cast<double>(ny);
    const double dz = lengthZ / static_cast<double>(nz);

    std::vector<Triangle3D> triangles;
    if (nx < 2 || ny < 2 || nz < 2) {
        return triangles;
    }

    for (std::size_t k = 0; k + 1 < nz; ++k) {
        for (std::size_t j = 0; j + 1 < ny; ++j) {
            for (std::size_t i = 0; i + 1 < nx; ++i) {
                // Cube corners 0..7, same convention as apps/unified_viewer's
                // cavity3d wireframe box: 0-3 the k face, 4-7 the k+1 face,
                // each quad ordered (i,j) (i+1,j) (i+1,j+1) (i,j+1).
                const std::array<Vector3, 8> c = {
                    cellCenter(i, j, k, dx, dy, dz),         cellCenter(i + 1, j, k, dx, dy, dz),
                    cellCenter(i + 1, j + 1, k, dx, dy, dz), cellCenter(i, j + 1, k, dx, dy, dz),
                    cellCenter(i, j, k + 1, dx, dy, dz),     cellCenter(i + 1, j, k + 1, dx, dy, dz),
                    cellCenter(i + 1, j + 1, k + 1, dx, dy, dz), cellCenter(i, j + 1, k + 1, dx, dy, dz),
                };
                const std::array<double, 8> f = {
                    field[i + j * nx + k * nx * ny],             field[(i + 1) + j * nx + k * nx * ny],
                    field[(i + 1) + (j + 1) * nx + k * nx * ny], field[i + (j + 1) * nx + k * nx * ny],
                    field[i + j * nx + (k + 1) * nx * ny],       field[(i + 1) + j * nx + (k + 1) * nx * ny],
                    field[(i + 1) + (j + 1) * nx + (k + 1) * nx * ny],
                    field[i + (j + 1) * nx + (k + 1) * nx * ny],
                };

                // Freudenthal decomposition into 6 tetrahedra sharing the
                // main diagonal c0-c6 -- see the header comment.
                static constexpr int kTets[6][4] = {
                    {0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6}, {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6},
                };
                for (const auto& tet : kTets) {
                    const std::array<Vector3, 4> p = {c[tet[0]], c[tet[1]], c[tet[2]], c[tet[3]]};
                    const std::array<double, 4> v = {f[tet[0]], f[tet[1]], f[tet[2]], f[tet[3]]};
                    processTetrahedron(p, v, isoValue, triangles);
                }
            }
        }
    }
    return triangles;
}

} // namespace aether::postprocessing
