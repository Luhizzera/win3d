#pragma once

#include "aether/mesh/TetrahedralMesh.hpp"
#include "aether/persistence/FieldArchive.hpp"

namespace aether::persistence {

// Module 11's mesh piece for the **unstructured** case, the one that had no
// persistence path at all.
//
// `GridArchive` covers StructuredGrid3D, which needs only six scalars
// because a uniform Cartesian grid is fully described by its corners and
// its counts. `engine/geometry`'s STL/OBJ I/O covers surface meshes. A
// tetrahedral volume mesh is neither: it carries real per-vertex positions
// and per-cell connectivity that nothing else in the project could write
// down, so a simulation on an imported geometry could be run but not
// *saved* -- the fields could be checkpointed while the mesh they were
// defined over could not, which makes the checkpoint unusable on its own.
//
// **Stored as connectivity, not as the generator's input.** The obvious
// alternative was to save the point cloud and re-tetrahedralize on load,
// which would be smaller. It is also wrong: the Delaunay result for a point
// set with co-spherical ties is not unique, so a reload could legitimately
// produce different cells, and the fields checkpointed alongside would then
// be indexed against a mesh that no longer matches them. Writing the cells
// explicitly is what makes the round trip an identity rather than a
// recomputation.
//
// Uses FieldArchive's array storage (a flat vertex array and a flat
// connectivity array) rather than inventing a second container format, so a
// mesh and the fields defined over it travel in one file and share its
// versioning.
void saveTetrahedralMesh(FieldArchive& archive, const mesh::TetrahedralMesh& mesh);

// Rebuilds the mesh a previous saveTetrahedralMesh() wrote.
//
// Throws std::runtime_error if the archive lacks the mesh arrays, or if the
// connectivity array's length is not a multiple of four -- a truncated file
// should fail where it is read rather than produce a mesh missing its last
// cell.
mesh::TetrahedralMesh loadTetrahedralMesh(const FieldArchive& archive);

} // namespace aether::persistence
