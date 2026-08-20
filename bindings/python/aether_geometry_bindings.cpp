#include "aether/geometry/ObjIO.hpp"
#include "aether/geometry/StlIO.hpp"
#include "aether/geometry/TriangleMesh.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::geometry;

PYBIND11_MODULE(aether_geometry_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's geometry layer (Module 2)";

    py::class_<TriangleMesh>(m, "TriangleMesh")
        .def(py::init<>())
        .def("add_vertex", &TriangleMesh::addVertex)
        .def("add_triangle", &TriangleMesh::addTriangle)
        .def("vertex_count", &TriangleMesh::vertexCount)
        .def("triangle_count", &TriangleMesh::triangleCount)
        .def("vertex", &TriangleMesh::vertex)
        // Exposes a triangle's three vertex indices as a plain tuple rather
        // than binding the Triangle struct itself -- this is the one piece
        // of TriangleMesh connectivity Python couldn't see at all before
        // (only vertex positions and counts were bound), and a tuple is all
        // a caller building facets for DelaunayTetrahedralization3D needs.
        .def("triangle",
             [](const TriangleMesh& mesh, std::size_t index) {
                 const Triangle& t = mesh.triangle(index);
                 return std::make_tuple(t.vertices[0], t.vertices[1], t.vertices[2]);
             })
        .def("face_normal", &TriangleMesh::faceNormal)
        .def("face_area", &TriangleMesh::faceArea)
        .def("surface_area", &TriangleMesh::surfaceArea)
        .def("volume", &TriangleMesh::volume)
        .def("weld_vertices", &TriangleMesh::weldVertices, py::arg("tolerance") = 1e-6)
        .def("remove_degenerate_triangles", &TriangleMesh::removeDegenerateTriangles,
             py::arg("area_tolerance") = 1e-12)
        .def("reorient_normals_outward", &TriangleMesh::reorientNormalsOutward)
        .def("find_boundary_edges",
             [](const TriangleMesh& mesh) {
                 std::vector<std::pair<std::size_t, std::size_t>> edges;
                 for (const auto& e : mesh.findBoundaryEdges()) {
                     edges.emplace_back(e.a, e.b);
                 }
                 return edges;
             })
        .def("is_watertight", &TriangleMesh::isWatertight);

    m.def("load_stl", &loadStl, py::arg("path"), py::arg("weld_tolerance") = 1e-6);
    m.def("save_stl_binary", &saveStlBinary, py::arg("mesh"), py::arg("path"));
    m.def("load_obj", &loadObj, py::arg("path"));
    m.def("save_obj", &saveObj, py::arg("mesh"), py::arg("path"));
}
