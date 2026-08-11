#pragma once

#include "aether/core/Mesh.hpp"

#include <cstddef>
#include <vector>

namespace aether::core {

// Cell-centered scalar field (e.g. pressure, temperature) defined over a
// Mesh. The Mesh is referenced, not owned or copied.
class ScalarField {
public:
    explicit ScalarField(const Mesh& mesh, double initialValue = 0.0);

    std::size_t size() const { return values_.size(); }

    double& operator[](std::size_t cellIndex) { return values_.at(cellIndex); }
    double operator[](std::size_t cellIndex) const { return values_.at(cellIndex); }

    void fill(double value);

    ScalarField operator+(const ScalarField& o) const;
    ScalarField operator-(const ScalarField& o) const;
    ScalarField operator*(double scalar) const;

private:
    const Mesh* mesh_;
    std::vector<double> values_;
};

} // namespace aether::core
