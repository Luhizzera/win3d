#pragma once

#include "aether/core/Vector3.hpp"

#include <cstddef>
#include <vector>

namespace aether::core {

// Minimal unstructured mesh data structure: vertices plus cells, each cell
// referencing the indices of its vertices. This is the data structure the
// physics/fields layer (ScalarField, VectorField) is defined over; mesh
// *generation* (tet/hex/prism, AMR) is a separate, later module.
class Mesh {
public:
    std::size_t addVertex(const Vector3& position);
    std::size_t addCell(std::vector<std::size_t> vertexIndices);

    std::size_t vertexCount() const { return vertices_.size(); }
    std::size_t cellCount() const { return cells_.size(); }

    const Vector3& vertex(std::size_t index) const { return vertices_.at(index); }
    const std::vector<std::size_t>& cell(std::size_t index) const { return cells_.at(index); }

    // Average of the cell's vertex positions. A placeholder for the volume-
    // weighted centroid formulas that a real mesh generator will provide.
    Vector3 cellCentroid(std::size_t cellIndex) const;

private:
    std::vector<Vector3> vertices_;
    std::vector<std::vector<std::size_t>> cells_;
};

} // namespace aether::core
