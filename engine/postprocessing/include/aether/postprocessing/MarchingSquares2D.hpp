#pragma once

#include "aether/core/Vector3.hpp"

#include <cstddef>
#include <vector>

namespace aether::postprocessing {

struct Segment2D {
    core::Vector3 a;
    core::Vector3 b;
};

// Module 7: iso-contour extraction from a 2D scalar field (marching
// squares) -- the classic first step toward marching cubes (3D iso-
// surfaces), deliberately scoped to 2D first per this project's usual
// incremental practice.
//
// field is an nx*ny array, cell-centered (same convention as every
// scalar field in this project: field(i,j) at ((i+0.5)*dx, (j+0.5)*dy)).
// For each of the (nx-1)*(ny-1) quads formed by 4 adjacent cell centers,
// classifies which corners are above isoValue, linearly interpolates the
// crossing point along each edge whose endpoints disagree, and connects
// them into 0, 1, or 2 line segments per quad -- derived directly from
// the crossed-edge count (always 0, 2, or 4 for a 4-cycle boundary) rather
// than a hardcoded 16-row case table, to avoid a transcription-error risk
// for a table that would otherwise have to be recalled from memory. The
// ambiguous 4-crossed-edge ("saddle") case is resolved with the average-
// corner-value heuristic (a standard, simple disambiguation -- not the
// full asymptotic decider some implementations use).
std::vector<Segment2D> marchingSquares2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY,
                                          const std::vector<double>& field, double isoValue);

} // namespace aether::postprocessing
