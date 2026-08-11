#include "aether/core/Mesh.hpp"

#include <stdexcept>

namespace aether::core {

std::size_t Mesh::addVertex(const Vector3& position) {
    vertices_.push_back(position);
    return vertices_.size() - 1;
}

std::size_t Mesh::addCell(std::vector<std::size_t> vertexIndices) {
    if (vertexIndices.empty()) {
        throw std::invalid_argument("Mesh::addCell: cell must reference at least one vertex");
    }
    for (std::size_t index : vertexIndices) {
        if (index >= vertices_.size()) {
            throw std::out_of_range("Mesh::addCell: vertex index out of range");
        }
    }
    cells_.push_back(std::move(vertexIndices));
    return cells_.size() - 1;
}

Vector3 Mesh::cellCentroid(std::size_t cellIndex) const {
    const auto& indices = cell(cellIndex);
    Vector3 sum;
    for (std::size_t index : indices) {
        sum += vertices_[index];
    }
    return sum / static_cast<double>(indices.size());
}

} // namespace aether::core
