#include "aether/core/ScalarField.hpp"

#include <stdexcept>

namespace aether::core {

ScalarField::ScalarField(const Mesh& mesh, double initialValue)
    : mesh_(&mesh), values_(mesh.cellCount(), initialValue) {}

void ScalarField::fill(double value) {
    for (double& v : values_) {
        v = value;
    }
}

ScalarField ScalarField::operator+(const ScalarField& o) const {
    if (size() != o.size()) {
        throw std::invalid_argument("ScalarField::operator+: size mismatch");
    }
    ScalarField result(*mesh_);
    for (std::size_t i = 0; i < size(); ++i) {
        result.values_[i] = values_[i] + o.values_[i];
    }
    return result;
}

ScalarField ScalarField::operator-(const ScalarField& o) const {
    if (size() != o.size()) {
        throw std::invalid_argument("ScalarField::operator-: size mismatch");
    }
    ScalarField result(*mesh_);
    for (std::size_t i = 0; i < size(); ++i) {
        result.values_[i] = values_[i] - o.values_[i];
    }
    return result;
}

ScalarField ScalarField::operator*(double scalar) const {
    ScalarField result(*mesh_);
    for (std::size_t i = 0; i < size(); ++i) {
        result.values_[i] = values_[i] * scalar;
    }
    return result;
}

} // namespace aether::core
