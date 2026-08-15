#include "aether/mesh/StructuredGrid3D.hpp"

#include <stdexcept>

namespace aether::mesh {

using aether::core::Vector3;

StructuredGrid3D::StructuredGrid3D(const Vector3& min, const Vector3& max, std::size_t nx,
                                    std::size_t ny, std::size_t nz)
    : min_(min), max_(max), nx_(nx), ny_(ny), nz_(nz) {
    if (nx == 0 || ny == 0 || nz == 0) {
        throw std::invalid_argument("StructuredGrid3D: nx, ny, nz must all be > 0");
    }
    spacing_ = Vector3((max.x - min.x) / static_cast<double>(nx), (max.y - min.y) / static_cast<double>(ny),
                       (max.z - min.z) / static_cast<double>(nz));
}

Vector3 StructuredGrid3D::cellCenter(std::size_t i, std::size_t j, std::size_t k) const {
    return Vector3(min_.x + (static_cast<double>(i) + 0.5) * spacing_.x,
                   min_.y + (static_cast<double>(j) + 0.5) * spacing_.y,
                   min_.z + (static_cast<double>(k) + 0.5) * spacing_.z);
}

bool StructuredGrid3D::hasNeighbor(std::size_t i, std::size_t j, std::size_t k, int di, int dj,
                                    int dk) const {
    const long long ni = static_cast<long long>(i) + di;
    const long long nj = static_cast<long long>(j) + dj;
    const long long nk = static_cast<long long>(k) + dk;
    return ni >= 0 && ni < static_cast<long long>(nx_) && nj >= 0 && nj < static_cast<long long>(ny_) &&
           nk >= 0 && nk < static_cast<long long>(nz_);
}

} // namespace aether::mesh
