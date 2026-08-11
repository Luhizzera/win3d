#include "aether/mesh/DelaunayTetrahedralization3D.hpp"
#include "aether/mesh/DelaunayTriangulation2D.hpp"
#include "aether/mesh/PolygonTriangulation2D.hpp"
#include "aether/mesh/StructuredGrid3D.hpp"

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
}
