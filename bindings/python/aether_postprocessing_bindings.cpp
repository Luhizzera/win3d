#include "aether/postprocessing/MarchingCubes3D.hpp"
#include "aether/postprocessing/MarchingSquares2D.hpp"
#include "aether/postprocessing/Streamline2D.hpp"
#include "aether/postprocessing/VtkWriter.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <map>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace aether::postprocessing;

PYBIND11_MODULE(aether_postprocessing_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's post-processing layer (Module 7)";

    py::class_<Streamline2D>(m, "Streamline2D")
        .def(py::init<std::size_t, std::size_t, double, double, const std::vector<double>&,
                      const std::vector<double>&, bool>(),
             py::arg("nx"), py::arg("ny"), py::arg("length_x"), py::arg("length_y"), py::arg("u"),
             py::arg("v"), py::arg("periodic"))
        .def("velocity_at", &Streamline2D::velocityAt, py::arg("x"), py::arg("y"))
        .def("trace", &Streamline2D::trace, py::arg("x0"), py::arg("y0"), py::arg("step_size"),
             py::arg("max_steps"));

    py::class_<Segment2D>(m, "Segment2D")
        .def_readonly("a", &Segment2D::a)
        .def_readonly("b", &Segment2D::b);

    m.def("marching_squares_2d", &marchingSquares2D, py::arg("nx"), py::arg("ny"), py::arg("length_x"),
          py::arg("length_y"), py::arg("field"), py::arg("iso_value"));

    py::class_<Triangle3D>(m, "Triangle3D")
        .def_readonly("a", &Triangle3D::a)
        .def_readonly("b", &Triangle3D::b)
        .def_readonly("c", &Triangle3D::c);

    m.def("marching_cubes_3d", &marchingCubes3D, py::arg("nx"), py::arg("ny"), py::arg("nz"),
          py::arg("length_x"), py::arg("length_y"), py::arg("length_z"), py::arg("field"),
          py::arg("iso_value"));

    // Exposed as plain dicts rather than by binding CellScalarField and
    // CellVectorField as classes: from Python the natural way to say "these
    // are my fields" is {"pressure": [...], "velocity": [...]}, and a caller
    // should not have to construct wrapper objects to hand over two lists.
    // The C++ structs stay the interface on that side, where naming a type
    // is what keeps a long argument list readable.
    m.def(
        "write_tetrahedral_mesh_vtk",
        [](const std::string& path, const aether::mesh::TetrahedralMesh& mesh,
           const std::map<std::string, std::vector<double>>& scalars,
           const std::map<std::string, std::vector<aether::core::Vector3>>& vectors) {
            std::vector<CellScalarField> scalarFields;
            scalarFields.reserve(scalars.size());
            for (const auto& [name, values] : scalars) {
                scalarFields.push_back(CellScalarField{name, values});
            }
            std::vector<CellVectorField> vectorFields;
            vectorFields.reserve(vectors.size());
            for (const auto& [name, values] : vectors) {
                vectorFields.push_back(CellVectorField{name, values});
            }
            writeTetrahedralMeshVtk(path, mesh, scalarFields, vectorFields);
        },
        py::arg("path"), py::arg("mesh"), py::arg("scalars") = std::map<std::string, std::vector<double>>{},
        py::arg("vectors") = std::map<std::string, std::vector<aether::core::Vector3>>{});
}
