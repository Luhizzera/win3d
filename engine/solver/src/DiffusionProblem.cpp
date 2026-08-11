#include "aether/solver/DiffusionProblem.hpp"

namespace aether::solver {

using aether::core::Vector3;
using aether::mesh::StructuredGrid3D;

DiffusionProblem::DiffusionProblem(const StructuredGrid3D& grid)
    : grid_(&grid), field_(grid.cellCount(), 0.0), isFixed_(grid.cellCount(), false) {}

void DiffusionProblem::setBoundaryValue(Face face, double value) {
    const std::size_t nx = grid_->nx();
    const std::size_t ny = grid_->ny();
    const std::size_t nz = grid_->nz();

    auto fix = [this, value](std::size_t idx) {
        field_[idx] = value;
        isFixed_[idx] = true;
    };

    switch (face) {
    case Face::XMin:
        for (std::size_t k = 0; k < nz; ++k) {
            for (std::size_t j = 0; j < ny; ++j) {
                fix(grid_->cellIndex(0, j, k));
            }
        }
        break;
    case Face::XMax:
        for (std::size_t k = 0; k < nz; ++k) {
            for (std::size_t j = 0; j < ny; ++j) {
                fix(grid_->cellIndex(nx - 1, j, k));
            }
        }
        break;
    case Face::YMin:
        for (std::size_t k = 0; k < nz; ++k) {
            for (std::size_t i = 0; i < nx; ++i) {
                fix(grid_->cellIndex(i, 0, k));
            }
        }
        break;
    case Face::YMax:
        for (std::size_t k = 0; k < nz; ++k) {
            for (std::size_t i = 0; i < nx; ++i) {
                fix(grid_->cellIndex(i, ny - 1, k));
            }
        }
        break;
    case Face::ZMin:
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                fix(grid_->cellIndex(i, j, 0));
            }
        }
        break;
    case Face::ZMax:
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                fix(grid_->cellIndex(i, j, nz - 1));
            }
        }
        break;
    }
}

void DiffusionProblem::setSourceTerm(double source) { source_ = source; }

std::array<std::pair<std::array<int, 3>, double>, 6> DiffusionProblem::neighborDirections() const {
    const Vector3 h = grid_->spacing();
    const double axWeight = 1.0 / (h.x * h.x);
    const double ayWeight = 1.0 / (h.y * h.y);
    const double azWeight = 1.0 / (h.z * h.z);

    return {{
        {{-1, 0, 0}, axWeight},
        {{1, 0, 0}, axWeight},
        {{0, -1, 0}, ayWeight},
        {{0, 1, 0}, ayWeight},
        {{0, 0, -1}, azWeight},
        {{0, 0, 1}, azWeight},
    }};
}

double DiffusionProblem::value(std::size_t i, std::size_t j, std::size_t k) const {
    return field_[grid_->cellIndex(i, j, k)];
}

} // namespace aether::solver
