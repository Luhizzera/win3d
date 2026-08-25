#include "aether/persistence/TetrahedralMeshArchive.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace aether::persistence {

namespace {

// Field names, in one place rather than repeated as literals at four call
// sites -- a typo in one of them would produce an archive that saves and
// then fails to load, with nothing pointing at which spelling was wrong.
constexpr const char* kVerticesField = "tetrahedral_mesh_vertices";
constexpr const char* kCellsField = "tetrahedral_mesh_cells";

} // namespace

void saveTetrahedralMesh(FieldArchive& archive, const mesh::TetrahedralMesh& mesh) {
    std::vector<double> vertices;
    vertices.reserve(mesh.vertexCount() * 3);
    for (std::size_t v = 0; v < mesh.vertexCount(); ++v) {
        const core::Vector3& p = mesh.vertex(v);
        vertices.push_back(p.x);
        vertices.push_back(p.y);
        vertices.push_back(p.z);
    }

    // FieldArchive stores doubles, so the connectivity indices ride as
    // doubles. That is exact rather than lossy for any mesh this project
    // can hold in memory: a double represents every integer up to 2^53
    // exactly, and a mesh with more than nine quadrillion vertices is not
    // the failure mode to design around.
    std::vector<double> cells;
    cells.reserve(mesh.cellCount() * 4);
    for (std::size_t c = 0; c < mesh.cellCount(); ++c) {
        for (const std::size_t index : mesh.cellVertices(c)) {
            cells.push_back(static_cast<double>(index));
        }
    }

    archive.setField(kVerticesField, std::move(vertices));
    archive.setField(kCellsField, std::move(cells));
}

mesh::TetrahedralMesh loadTetrahedralMesh(const FieldArchive& archive) {
    if (!archive.hasField(kVerticesField) || !archive.hasField(kCellsField)) {
        throw std::runtime_error(
            "loadTetrahedralMesh: the archive carries no tetrahedral mesh (saved by "
            "saveTetrahedralMesh?)");
    }
    const std::vector<double>& vertexData = archive.field(kVerticesField);
    const std::vector<double>& cellData = archive.field(kCellsField);

    if (vertexData.size() % 3 != 0) {
        throw std::runtime_error(
            "loadTetrahedralMesh: the vertex array's length is not a multiple of 3");
    }
    if (cellData.size() % 4 != 0) {
        throw std::runtime_error(
            "loadTetrahedralMesh: the connectivity array's length is not a multiple of 4");
    }

    std::vector<core::Vector3> vertices;
    vertices.reserve(vertexData.size() / 3);
    for (std::size_t i = 0; i + 2 < vertexData.size(); i += 3) {
        vertices.push_back(core::Vector3(vertexData[i], vertexData[i + 1], vertexData[i + 2]));
    }

    std::vector<std::array<std::size_t, 4>> cells;
    cells.reserve(cellData.size() / 4);
    for (std::size_t i = 0; i + 3 < cellData.size(); i += 4) {
        std::array<std::size_t, 4> cell{};
        for (std::size_t k = 0; k < 4; ++k) {
            const double value = cellData[i + k];
            // A negative or fractional index means the array is not
            // connectivity at all -- caught here rather than wrapping into
            // a huge size_t and being reported as an out-of-range vertex.
            if (value < 0.0 || value != std::floor(value)) {
                throw std::runtime_error(
                    "loadTetrahedralMesh: the connectivity array holds a value that is not a "
                    "non-negative whole number");
            }
            cell[k] = static_cast<std::size_t>(value);
        }
        cells.push_back(cell);
    }

    // fromCells does the bounds checking against the vertex list, and
    // rebuilds face connectivity, volumes and centroids exactly as the
    // generator path does.
    return mesh::TetrahedralMesh::fromCells(vertices, cells);
}

} // namespace aether::persistence
