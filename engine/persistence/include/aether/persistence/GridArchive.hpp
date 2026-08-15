#pragma once

#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/persistence/FieldArchive.hpp"

namespace aether::persistence {

// Module 11's "mesh" piece for the structured grid case: StructuredGrid3D
// carries no per-cell data of its own (just min/max corner + nx/ny/nz --
// see its own header), so unlike a solver's checkpoint this doesn't need
// FieldArchive's array storage at all, only its scalar metadata map. These
// two free functions exist so callers don't have to remember the six
// metadata key names by hand every time a grid needs to travel alongside a
// solver's field data in the same archive.
//
// (TriangleMesh, the other mesh representation in this project, already
// has its own dedicated STL/OBJ import/export in engine/geometry -- that
// already covers Module 11's "mesh" persistence need for surface meshes;
// this is only for the structured volumetric grid, which had no
// persistence path at all before.)
//
// saveGrid() and loadGrid() round-trip via the exact same constructor
// arguments (min corner, max corner, nx, ny, nz) StructuredGrid3D's own
// constructor takes, rather than trying to reconstruct max from the
// stored spacing -- spacing is derived from max via division, and
// re-deriving max from spacing via multiplication is not guaranteed to
// invert that division bit-for-bit, so storing max directly is the exact,
// not approximate, choice.
void saveGrid(FieldArchive& archive, const mesh::StructuredGrid3D& grid);
mesh::StructuredGrid3D loadGrid(const FieldArchive& archive);

} // namespace aether::persistence
