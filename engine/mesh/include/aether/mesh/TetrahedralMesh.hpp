#pragma once

#include "aether/core/Mesh.hpp"
#include "aether/core/Vector3.hpp"
#include "aether/mesh/DelaunayTetrahedralization3D.hpp"

#include <array>
#include <cstddef>
#include <vector>

namespace aether::mesh {

// Face connectivity over a tetrahedral mesh -- the piece that was missing
// between this project's mesh *generation* and any finite-volume solver.
//
// **Why this class had to exist before unstructured FVM could.**
// DelaunayTetrahedralization3D produces tetrahedra as vertex quadruples,
// and core::Mesh stores cells as vertex lists. Neither knows which cells
// share a face, and a finite-volume method is defined entirely in terms of
// fluxes *through faces*: it needs, for every face, the two cells on
// either side, the face's area and normal, and the vector between the two
// cell centroids. Verified directly in the code before starting this:
// nothing in engine/solver referenced the tetrahedralizer at all, and all
// five solvers took only StructuredGrid3D -- so the sophisticated
// unstructured meshing built in Module 3 had never fed a solver.
//
// **What is deliberately exact here.** A tetrahedron's centroid is exactly
// the average of its four vertices (unlike a general polyhedron, where the
// vertex average is only an approximation), and its volume is exactly
// |(b-a) . ((c-a) x (d-a))| / 6. So no geometric quantity in this class is
// an approximation -- which is what makes the invariants below checkable
// as identities rather than as tolerances.
class TetrahedralMesh {
public:
    static constexpr std::size_t kNoNeighbour = static_cast<std::size_t>(-1);

    struct Face {
        std::array<std::size_t, 3> vertices;
        std::size_t owner;     // the cell the area vector points away from
        std::size_t neighbour; // kNoNeighbour on a boundary face
        // Outward normal scaled by the face's area, directed owner ->
        // neighbour (and out of the domain on a boundary face). Storing the
        // product rather than normal and area separately is what makes the
        // divergence-theorem identity below exact: it avoids normalizing
        // and re-scaling, which would introduce roundoff.
        core::Vector3 areaVector;
        core::Vector3 centroid;
    };

    // Builds the connectivity of an already-computed tetrahedralization.
    // Tetrahedra with (numerically) zero volume are skipped rather than
    // producing degenerate faces -- they carry no flux and would divide by
    // zero downstream.
    static TetrahedralMesh fromTetrahedralization(const DelaunayTetrahedralization3D& tetrahedralization);

    // The same construction from plain vertex positions and cell
    // connectivity, for a mesh that did not come from this project's own
    // Delaunay generator.
    //
    // **This existed as a gap, not a choice.** Until it was added, the only
    // route to a TetrahedralMesh was through DelaunayTetrahedralization3D,
    // which meant a mesh could not be *reloaded* from a checkpoint (see
    // engine/persistence/TetrahedralMeshArchive) or imported from any other
    // mesher without re-running a tetrahedralization that would not
    // necessarily reproduce the same cells. fromTetrahedralization() is now
    // a thin adapter over this.
    //
    // The winding of each cell is checked and corrected internally, exactly
    // as it already was: a tetrahedron from the Delaunay generator is
    // positively oriented by construction, and one from a file is whatever
    // the file says.
    //
    // Throws std::invalid_argument if a cell references a vertex index that
    // does not exist -- a corrupt file should fail at the boundary, not
    // read out of bounds somewhere downstream.
    static TetrahedralMesh fromCells(const std::vector<core::Vector3>& vertices,
                                      const std::vector<std::array<std::size_t, 4>>& cells);

    // **Vertices and cells live in a core::Mesh, not in a second copy here.**
    // That was DIVIDA_TECNICA.md 2.2: this class was built in Fase 2.1
    // alongside core::Mesh rather than on top of it, so the engine carried two
    // unstructured mesh representations that did not know about each other --
    // and since ScalarField/VectorField are defined over core::Mesh, every
    // unstructured solver here had to carry raw std::vector<double> instead of
    // the engine's own field types. Two representations of one thing is the
    // situation item 2.1 dealt with for solvers, one level down.
    //
    // core::Mesh is the canonical one because it is what the field layer
    // already speaks. What this class adds on top is *connectivity* -- which
    // cells share a face, and that face's geometry -- which is the part a
    // finite-volume method needs and a vertex-and-cell list cannot express.
    //
    // The two agree **exactly**, not approximately, and that is worth stating
    // because it is not true in general: core::Mesh::cellCentroid() is the
    // average of a cell's vertices, and for a tetrahedron that *is* the
    // centroid. For any other cell shape it would be an approximation and this
    // identity would need re-examining.
    const core::Mesh& coreMesh() const { return mesh_; }

    std::size_t cellCount() const { return cellVolumes_.size(); }
    std::size_t faceCount() const { return faces_.size(); }
    std::size_t vertexCount() const { return mesh_.vertexCount(); }

    const Face& face(std::size_t index) const { return faces_.at(index); }
    const core::Vector3& vertex(std::size_t index) const { return mesh_.vertex(index); }
    // A tetrahedron's four vertices, as core::Mesh stores them.
    const std::vector<std::size_t>& cellVertices(std::size_t cell) const { return mesh_.cell(cell); }
    double cellVolume(std::size_t cell) const { return cellVolumes_.at(cell); }
    const core::Vector3& cellCentroid(std::size_t cell) const { return cellCentroids_.at(cell); }
    const std::vector<std::size_t>& cellFaces(std::size_t cell) const { return cellFaces_.at(cell); }

    bool isBoundaryFace(std::size_t index) const { return faces_.at(index).neighbour == kNoNeighbour; }
    std::size_t boundaryFaceCount() const;

    // The area vector of `face` as seen from `cell`, i.e. flipped when
    // `cell` is the neighbour rather than the owner. Solvers want the
    // outward-from-me direction; storing one canonical direction and
    // flipping on request keeps the two sides exactly antisymmetric.
    core::Vector3 outwardAreaVector(std::size_t cell, std::size_t faceIndex) const;

    // Sum of a cell's outward area vectors. For any closed polyhedron this
    // is identically zero -- the discrete divergence theorem. Returned
    // rather than asserted so tests can measure how close to zero it lands
    // and so a caller can check a mesh it did not build. A nonzero result
    // means a normal is misoriented or a face is missing, which are exactly
    // the two bugs this class could plausibly have.
    core::Vector3 cellAreaVectorSum(std::size_t cell) const;

    double totalVolume() const;

private:
    core::Mesh mesh_;
    std::vector<Face> faces_;
    std::vector<double> cellVolumes_;
    std::vector<core::Vector3> cellCentroids_;
    std::vector<std::vector<std::size_t>> cellFaces_;
};

} // namespace aether::mesh
