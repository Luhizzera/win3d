#include "aether/persistence/GridArchive.hpp"

namespace aether::persistence {

using aether::core::Vector3;
using aether::mesh::StructuredGrid3D;

void saveGrid(FieldArchive& archive, const StructuredGrid3D& grid) {
    archive.setMetadata("grid_min_x", grid.min().x);
    archive.setMetadata("grid_min_y", grid.min().y);
    archive.setMetadata("grid_min_z", grid.min().z);
    archive.setMetadata("grid_max_x", grid.max().x);
    archive.setMetadata("grid_max_y", grid.max().y);
    archive.setMetadata("grid_max_z", grid.max().z);
    archive.setMetadata("grid_nx", static_cast<double>(grid.nx()));
    archive.setMetadata("grid_ny", static_cast<double>(grid.ny()));
    archive.setMetadata("grid_nz", static_cast<double>(grid.nz()));
}

StructuredGrid3D loadGrid(const FieldArchive& archive) {
    const Vector3 min(archive.metadata("grid_min_x"), archive.metadata("grid_min_y"), archive.metadata("grid_min_z"));
    const Vector3 max(archive.metadata("grid_max_x"), archive.metadata("grid_max_y"), archive.metadata("grid_max_z"));
    const auto nx = static_cast<std::size_t>(archive.metadata("grid_nx"));
    const auto ny = static_cast<std::size_t>(archive.metadata("grid_ny"));
    const auto nz = static_cast<std::size_t>(archive.metadata("grid_nz"));
    return StructuredGrid3D(min, max, nx, ny, nz);
}

} // namespace aether::persistence
