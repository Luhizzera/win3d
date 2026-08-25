#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"

#include <string>
#include <vector>

namespace aether::postprocessing {

// Module 7: writes a TetrahedralMesh and the fields defined over it to a
// VTK legacy ASCII unstructured-grid file (.vtk).
//
// **Why an exporter is the right first answer to "how do I look at an
// unstructured result", ahead of building a viewer for it.** Every other
// solver in this project has a built-in viewer mode, and those were the
// right call: a structured grid's field is a dense array, trivially
// rendered as a heatmap or an arrow lattice. An unstructured tetrahedral
// result is a different problem -- slicing, thresholding, streamline
// seeding and iso-surfacing over an irregular mesh are exactly the
// operations a general post-processor exists to provide, and ParaView and
// VisIt both read this format directly. Writing it costs a few hundred
// lines here against re-implementing a post-processor's worth of
// interaction inside a Win32/GL demo app.
//
// **The legacy ASCII format, not XML .vtu**, deliberately: it is a stable,
// fully-specified, human-readable text layout that needs no XML writer, no
// base64, no compression library and no schema version negotiation, and
// every tool that reads VTK at all reads it. The cost is file size and
// parse speed, which matter for a production pipeline dumping thousands of
// timesteps and do not matter for looking at a result.
//
// Fields are **cell-centred**, matching where this engine's finite-volume
// unknowns actually live: UnstructuredCavitySolver3D stores velocity and
// pressure per cell, and writing them as point data would mean inventing
// an interpolation to the vertices and then presenting the interpolated
// values as if they were the solution.
struct CellScalarField {
    std::string name;
    std::vector<double> values; // one per cell
};

struct CellVectorField {
    std::string name;
    std::vector<core::Vector3> values; // one per cell
};

// Writes `mesh` plus any cell fields to `path`.
//
// Throws std::invalid_argument if a field's length does not match the
// mesh's cell count, or if a field name is empty or contains whitespace --
// the format separates tokens by whitespace, so a name with a space in it
// produces a file that parses as something else entirely rather than
// failing loudly. Throws std::runtime_error if the file cannot be opened.
//
// Coordinates and field values are written with enough decimal digits to
// round-trip a double exactly (17 significant digits), because a result
// that changes when it is written out is not the result.
void writeTetrahedralMeshVtk(const std::string& path, const mesh::TetrahedralMesh& mesh,
                              const std::vector<CellScalarField>& scalars = {},
                              const std::vector<CellVectorField>& vectors = {});

} // namespace aether::postprocessing
