#include "aether/mesh/DelaunayTetrahedralization3D.hpp"
#include "aether/mesh/DelaunayTriangulation2D.hpp"
#include "aether/mesh/PolygonTriangulation2D.hpp"
#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::mesh;

PYBIND11_MODULE(aether_mesh_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's mesh layer (Module 3)";

    py::class_<StructuredGrid3D>(m, "StructuredGrid3D")
        .def(py::init<const aether::core::Vector3&, const aether::core::Vector3&, std::size_t, std::size_t,
                      std::size_t>(),
             py::arg("min"), py::arg("max"), py::arg("nx"), py::arg("ny"), py::arg("nz"))
        .def("nx", &StructuredGrid3D::nx)
        .def("ny", &StructuredGrid3D::ny)
        .def("nz", &StructuredGrid3D::nz)
        .def("cell_count", &StructuredGrid3D::cellCount)
        .def("min", &StructuredGrid3D::min)
        .def("max", &StructuredGrid3D::max)
        .def("spacing", &StructuredGrid3D::spacing)
        .def("cell_volume", &StructuredGrid3D::cellVolume)
        .def("cell_index", &StructuredGrid3D::cellIndex)
        .def("cell_center", &StructuredGrid3D::cellCenter)
        .def("has_neighbor", &StructuredGrid3D::hasNeighbor);

    py::class_<DelaunayTriangulation2D>(m, "DelaunayTriangulation2D")
        .def(py::init<>())
        .def("add_point", &DelaunayTriangulation2D::addPoint, py::arg("x"), py::arg("y"))
        .def("triangulate", &DelaunayTriangulation2D::triangulate)
        .def("point_count", &DelaunayTriangulation2D::pointCount)
        .def("triangle_count", &DelaunayTriangulation2D::triangleCount)
        .def("point", &DelaunayTriangulation2D::point)
        .def("triangle", [](const DelaunayTriangulation2D& tri, std::size_t i) {
            return tri.triangle(i).vertices;
        })
        .def("satisfies_delaunay_property", &DelaunayTriangulation2D::satisfiesDelaunayProperty,
             py::arg("tolerance") = 1e-7);

    py::class_<PolygonTriangulation2D>(m, "PolygonTriangulation2D")
        .def(py::init<>())
        .def("add_vertex", &PolygonTriangulation2D::addVertex, py::arg("x"), py::arg("y"))
        .def("triangulate", &PolygonTriangulation2D::triangulate)
        .def("vertex_count", &PolygonTriangulation2D::vertexCount)
        .def("triangle_count", &PolygonTriangulation2D::triangleCount)
        .def("vertex", &PolygonTriangulation2D::vertex)
        .def("triangle", [](const PolygonTriangulation2D& tri, std::size_t i) {
            return tri.triangle(i).vertices;
        })
        .def("polygon_area", &PolygonTriangulation2D::polygonArea)
        .def("is_locally_delaunay", &PolygonTriangulation2D::isLocallyDelaunay, py::arg("tolerance") = 1e-7);

    py::class_<DelaunayTetrahedralization3D::FacetRecoveryResult>(m, "FacetRecoveryResult")
        .def_readonly("recovered_facets", &DelaunayTetrahedralization3D::FacetRecoveryResult::recoveredFacets)
        .def_readonly("unrecovered", &DelaunayTetrahedralization3D::FacetRecoveryResult::unrecovered);

    py::class_<DelaunayTetrahedralization3D>(m, "DelaunayTetrahedralization3D")
        .def(py::init<>())
        .def("add_point", &DelaunayTetrahedralization3D::addPoint, py::arg("x"), py::arg("y"), py::arg("z"))
        .def("tetrahedralize", &DelaunayTetrahedralization3D::tetrahedralize)
        .def("insert_steiner_point", &DelaunayTetrahedralization3D::insertSteinerPoint, py::arg("x"),
             py::arg("y"), py::arg("z"))
        .def("missing_facets", &DelaunayTetrahedralization3D::missingFacets, py::arg("facets"))
        .def("recover_facets", &DelaunayTetrahedralization3D::recoverFacets, py::arg("facets"),
             py::arg("max_rounds") = 4)
        .def("remove_region", &DelaunayTetrahedralization3D::removeRegion, py::arg("seed"), py::arg("walls"))
        .def("point_count", &DelaunayTetrahedralization3D::pointCount)
        .def("tetrahedron_count", &DelaunayTetrahedralization3D::tetrahedronCount)
        .def("point", &DelaunayTetrahedralization3D::point)
        .def("tetrahedron", [](const DelaunayTetrahedralization3D& tet, std::size_t i) {
            return tet.tetrahedron(i).vertices;
        })
        .def("satisfies_delaunay_property", &DelaunayTetrahedralization3D::satisfiesDelaunayProperty,
             py::arg("tolerance") = 1e-6);

    // Face connectivity, the piece the unstructured solvers consume. Bound
    // together with them (see aether_solver_bindings.cpp) because a solver
    // without a way to build its mesh from Python would expose the layer
    // without making it usable.
    py::class_<TetrahedralMesh::Face>(m, "TetrahedralFace")
        .def_readonly("vertices", &TetrahedralMesh::Face::vertices)
        .def_readonly("owner", &TetrahedralMesh::Face::owner)
        // kNoNeighbour (the largest size_t) on a boundary face -- compare
        // against TetrahedralMesh.NO_NEIGHBOUR rather than to a literal, or
        // just ask is_boundary_face().
        .def_readonly("neighbour", &TetrahedralMesh::Face::neighbour)
        .def_readonly("area_vector", &TetrahedralMesh::Face::areaVector)
        .def_readonly("centroid", &TetrahedralMesh::Face::centroid);

    py::class_<TetrahedralMesh>(m, "TetrahedralMesh")
        .def_static("from_tetrahedralization", &TetrahedralMesh::fromTetrahedralization,
                    py::arg("tetrahedralization"))
        .def_readonly_static("NO_NEIGHBOUR", &TetrahedralMesh::kNoNeighbour)
        .def("cell_count", &TetrahedralMesh::cellCount)
        .def("face_count", &TetrahedralMesh::faceCount)
        .def("vertex_count", &TetrahedralMesh::vertexCount)
        // Returned by reference into the mesh's own storage, so the mesh has
        // to outlive the Face handed back.
        .def("face", &TetrahedralMesh::face, py::arg("index"),
             py::return_value_policy::reference_internal)
        .def("vertex", &TetrahedralMesh::vertex, py::arg("index"))
        .def("cell_volume", &TetrahedralMesh::cellVolume, py::arg("cell"))
        .def("cell_centroid", &TetrahedralMesh::cellCentroid, py::arg("cell"))
        .def("cell_faces", &TetrahedralMesh::cellFaces, py::arg("cell"))
        .def("is_boundary_face", &TetrahedralMesh::isBoundaryFace, py::arg("index"))
        .def("boundary_face_count", &TetrahedralMesh::boundaryFaceCount)
        .def("outward_area_vector", &TetrahedralMesh::outwardAreaVector, py::arg("cell"),
             py::arg("face_index"))
        // Identically zero for any closed polyhedron -- the discrete
        // divergence theorem. Exposed because it is how a caller checks a
        // mesh it did not build: a nonzero result means a normal is
        // misoriented or a face is missing.
        .def("cell_area_vector_sum", &TetrahedralMesh::cellAreaVectorSum, py::arg("cell"))
        .def("total_volume", &TetrahedralMesh::totalVolume);
}
