#include "aether/postprocessing/VtkWriter.hpp"

#include <cctype>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace aether::postprocessing {

namespace {

// VTK's own cell-type code for a linear tetrahedron. Named rather than
// written as a bare 10 at the point of use, because a wrong number here
// produces a file that loads without complaint and displays nonsense.
constexpr int kVtkTetraCellType = 10;

void validateName(const std::string& name) {
    if (name.empty()) {
        throw std::invalid_argument("writeTetrahedralMeshVtk: field name must not be empty");
    }
    for (const char c : name) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            throw std::invalid_argument("writeTetrahedralMeshVtk: field name '" + name +
                                         "' contains whitespace, which the legacy VTK format uses "
                                         "as its token separator");
        }
    }
}

} // namespace

void writeTetrahedralMeshVtk(const std::string& path, const mesh::TetrahedralMesh& mesh,
                              const std::vector<CellScalarField>& scalars,
                              const std::vector<CellVectorField>& vectors) {
    const std::size_t cellCount = mesh.cellCount();
    for (const CellScalarField& field : scalars) {
        validateName(field.name);
        if (field.values.size() != cellCount) {
            throw std::invalid_argument("writeTetrahedralMeshVtk: scalar field '" + field.name +
                                         "' has a different length than the mesh's cell count");
        }
    }
    for (const CellVectorField& field : vectors) {
        validateName(field.name);
        if (field.values.size() != cellCount) {
            throw std::invalid_argument("writeTetrahedralMeshVtk: vector field '" + field.name +
                                         "' has a different length than the mesh's cell count");
        }
    }

    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("writeTetrahedralMeshVtk: cannot open '" + path + "' for writing");
    }
    // Enough digits to round-trip an IEEE754 double exactly.
    out.precision(std::numeric_limits<double>::max_digits10);

    out << "# vtk DataFile Version 3.0\n";
    out << "Aether CFD Engine unstructured result\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    const std::size_t vertexCount = mesh.vertexCount();
    out << "POINTS " << vertexCount << " double\n";
    for (std::size_t v = 0; v < vertexCount; ++v) {
        const core::Vector3& p = mesh.vertex(v);
        out << p.x << ' ' << p.y << ' ' << p.z << '\n';
    }

    // Each tetrahedron is written as its vertex count followed by its four
    // indices, so the total integer count the header declares is 5 per cell.
    out << "CELLS " << cellCount << ' ' << cellCount * 5 << '\n';
    for (std::size_t c = 0; c < cellCount; ++c) {
        const std::vector<std::size_t>& vertices = mesh.cellVertices(c);
        out << "4";
        for (const std::size_t index : vertices) {
            out << ' ' << index;
        }
        out << '\n';
    }

    out << "CELL_TYPES " << cellCount << '\n';
    for (std::size_t c = 0; c < cellCount; ++c) {
        out << kVtkTetraCellType << '\n';
    }

    if (!scalars.empty() || !vectors.empty()) {
        out << "CELL_DATA " << cellCount << '\n';
        for (const CellScalarField& field : scalars) {
            out << "SCALARS " << field.name << " double 1\n";
            out << "LOOKUP_TABLE default\n";
            for (const double value : field.values) {
                out << value << '\n';
            }
        }
        for (const CellVectorField& field : vectors) {
            out << "VECTORS " << field.name << " double\n";
            for (const core::Vector3& value : field.values) {
                out << value.x << ' ' << value.y << ' ' << value.z << '\n';
            }
        }
    }

    out.flush();
    if (!out) {
        throw std::runtime_error("writeTetrahedralMeshVtk: failed while writing '" + path + "'");
    }
}

} // namespace aether::postprocessing
