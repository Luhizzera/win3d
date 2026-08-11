#pragma once

#include "aether/mesh/StructuredGrid3D.hpp"

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace aether::solver {

// Shared state and discretization for anything built on a cell-centered
// scalar field phi over a StructuredGrid3D obeying nabla^2(phi) = -source
// (steady) or d(phi)/dt = nabla^2(phi) + source (transient) -- the
// boundary-condition bookkeeping, per-axis anisotropic stencil weights,
// and field storage are identical either way. SteadyDiffusionSolver and
// TransientDiffusionSolver each add their own time-stepping/solve method
// on top of this.
//
// Boundary conditions: setBoundaryValue() fixes an entire face's cell
// *layer* directly to a given value (not applied via a ghost-cell wall
// treatment at the true domain edge, so e.g. the exact discrete steady
// solution along an axis with both ends fixed and source == 0 is linear
// in cell *index*, not in physical position measured from the domain
// edges -- the two differ by half a cell width at each end). A face with
// no setBoundaryValue() call is insulated (zero-gradient), which falls out
// naturally of only averaging over neighbors that exist. Corner/edge cells
// shared by more than one face take whichever face was set last. A proper
// wall-distance treatment is future work.
//
// Grid spacing may differ per axis (hx, hy, hz need not be equal): each
// neighbor's contribution is weighted by 1/h_axis^2, the correct finite-
// volume discretization of the Laplacian on an anisotropic Cartesian grid.
class DiffusionProblem {
public:
    enum class Face { XMin, XMax, YMin, YMax, ZMin, ZMax };

    explicit DiffusionProblem(const mesh::StructuredGrid3D& grid);
    virtual ~DiffusionProblem() = default;

    void setBoundaryValue(Face face, double value);

    // Uniform volumetric source term (already divided by whatever
    // conductivity/viscosity applies -- see derived classes). Defaults to
    // 0.
    void setSourceTerm(double source);

    double value(std::size_t i, std::size_t j, std::size_t k) const;

protected:
    // The 6 face-neighbor directions paired with the per-axis weight
    // (1/h_axis^2) their contribution carries.
    std::array<std::pair<std::array<int, 3>, double>, 6> neighborDirections() const;

    const mesh::StructuredGrid3D* grid_;
    std::vector<double> field_;
    std::vector<bool> isFixed_;
    double source_ = 0.0;
};

} // namespace aether::solver
