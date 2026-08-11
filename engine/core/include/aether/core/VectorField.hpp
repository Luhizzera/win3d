#pragma once

#include "aether/core/Mesh.hpp"
#include "aether/core/Vector3.hpp"

#include <cstddef>
#include <vector>

namespace aether::core {

// Cell-centered vector field (e.g. velocity) defined over a Mesh. The Mesh
// is referenced, not owned or copied.
class VectorField {
public:
    explicit VectorField(const Mesh& mesh, const Vector3& initialValue = Vector3{});

    std::size_t size() const { return values_.size(); }

    Vector3& operator[](std::size_t cellIndex) { return values_.at(cellIndex); }
    const Vector3& operator[](std::size_t cellIndex) const { return values_.at(cellIndex); }

    void fill(const Vector3& value);

    VectorField operator+(const VectorField& o) const;
    VectorField operator-(const VectorField& o) const;
    VectorField operator*(double scalar) const;

private:
    const Mesh* mesh_;
    std::vector<Vector3> values_;
};

} // namespace aether::core
