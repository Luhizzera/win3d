#include "aether/core/VectorField.hpp"

#include <stdexcept>

namespace aether::core {

VectorField::VectorField(const Mesh& mesh, const Vector3& initialValue)
    : mesh_(&mesh), values_(mesh.cellCount(), initialValue) {}

void VectorField::fill(const Vector3& value) {
    for (Vector3& v : values_) {
        v = value;
    }
}

VectorField VectorField::operator+(const VectorField& o) const {
    if (size() != o.size()) {
        throw std::invalid_argument("VectorField::operator+: size mismatch");
    }
    VectorField result(*mesh_);
    for (std::size_t i = 0; i < size(); ++i) {
        result.values_[i] = values_[i] + o.values_[i];
    }
    return result;
}

VectorField VectorField::operator-(const VectorField& o) const {
    if (size() != o.size()) {
        throw std::invalid_argument("VectorField::operator-: size mismatch");
    }
    VectorField result(*mesh_);
    for (std::size_t i = 0; i < size(); ++i) {
        result.values_[i] = values_[i] - o.values_[i];
    }
    return result;
}

VectorField VectorField::operator*(double scalar) const {
    VectorField result(*mesh_);
    for (std::size_t i = 0; i < size(); ++i) {
        result.values_[i] = values_[i] * scalar;
    }
    return result;
}

} // namespace aether::core
