#include "aether/postprocessing/MarchingCubes3D.hpp"
#include "aether/postprocessing/MarchingSquares2D.hpp"
#include "aether/postprocessing/Streamline2D.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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
}
