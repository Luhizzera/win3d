#pragma once

#include "aether/core/Vector3.hpp"

#include <cstddef>
#include <vector>

namespace aether::postprocessing {

struct Triangle3D {
    core::Vector3 a;
    core::Vector3 b;
    core::Vector3 c;
};

// The 3D iso-surface extractor MarchingSquares2D's own header comment
// already flagged as "the classic next step" -- closes that thread, and is
// the prerequisite the cavity3d viewer mode's comment names for eventually
// rendering the 3D turbulence closures' scalar fields (nu_t, k,
// epsilon/omega) as surfaces instead of arrows.
//
// **Implemented via marching *tetrahedra*, not the classic 256-case cube
// table.** Recalling that table from memory (or the 15 canonical cases +
// rotations/reflections/complements it's usually compressed into) would be
// exactly the kind of transcription-error risk this project's practice
// avoids -- the same reasoning MarchingSquares2D's own comment gives for
// deriving its 4-crossing cases from an edge-crossing count instead of a
// 16-row table. Marching tetrahedra sidesteps the problem entirely: each
// cube is split into 6 tetrahedra sharing the main diagonal (the standard
// Freudenthal decomposition; validated in this class's test by summing the
// 6 sub-volumes of a unit cube and checking they total exactly 1.0, the
// same "measure an exact invariant" discipline
// DelaunayTetrahedralization3D's own bipyramid-volume test uses), and a
// tetrahedron has only 2^4 = 16 corner-sign combinations -- few enough, and
// simple enough (no ambiguous-face cases exist for a tetrahedron, unlike a
// cube face), to derive per-case from the count of "above-isoValue"
// corners rather than reciting anything:
//   0 or 4 corners above: no crossing, no triangle.
//   1 or 3 corners above: the lone differing corner's 3 edges to the other
//     three are exactly the crossed edges -- one triangle.
//   2 corners above, 2 below: the 4 edges connecting the "above" pair to
//     the "below" pair are exactly the crossed edges, forming a
//     quadrilateral (traced in cyclic order across the bipartite pairing)
//     split into 2 triangles.
//
// field is an nx*ny*nz array, cell-centered (same convention as every
// scalar field in this project: field(i,j,k) at
// ((i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz)). Triangles are wound with their
// normal pointing from below isoValue toward above isoValue (verified,
// not just asserted, by the class's own sphere test: generated normals
// point away from the sphere's center, matching a field that increases
// outward).
std::vector<Triangle3D> marchingCubes3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                         double lengthY, double lengthZ, const std::vector<double>& field,
                                         double isoValue);

} // namespace aether::postprocessing
