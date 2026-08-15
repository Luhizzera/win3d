#pragma once

#include "aether/core/Vector3.hpp"

#include <cstddef>

namespace aether::mesh {

// Uniform axis-aligned structured grid of hexahedral cells: nx * ny * nz
// cells filling the box [min, max]. Cell (i, j, k) satisfies 0 <= i < nx,
// 0 <= j < ny, 0 <= k < nz.
//
// This is a deliberately simpler, structured-only mesh representation used
// to get the first FVM solver (Module 5) working with O(1) neighbor
// lookups. It does not replace aether::core::Mesh: a general unstructured
// mesh generator (tet/hex/prism, AMR) producing that generic representation
// is separate, later work.
class StructuredGrid3D {
public:
    StructuredGrid3D(const core::Vector3& min, const core::Vector3& max, std::size_t nx,
                      std::size_t ny, std::size_t nz);

    std::size_t nx() const { return nx_; }
    std::size_t ny() const { return ny_; }
    std::size_t nz() const { return nz_; }
    std::size_t cellCount() const { return nx_ * ny_ * nz_; }

    const core::Vector3& min() const { return min_; }
    const core::Vector3& max() const { return max_; }
    const core::Vector3& spacing() const { return spacing_; }
    double cellVolume() const { return spacing_.x * spacing_.y * spacing_.z; }

    std::size_t cellIndex(std::size_t i, std::size_t j, std::size_t k) const {
        return i + j * nx_ + k * nx_ * ny_;
    }

    core::Vector3 cellCenter(std::size_t i, std::size_t j, std::size_t k) const;

    // True if (i + di, j + dj, k + dk) is inside the grid.
    bool hasNeighbor(std::size_t i, std::size_t j, std::size_t k, int di, int dj, int dk) const;

private:
    core::Vector3 min_;
    core::Vector3 max_;
    core::Vector3 spacing_;
    std::size_t nx_;
    std::size_t ny_;
    std::size_t nz_;
};

} // namespace aether::mesh
