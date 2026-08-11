#include "aether/core/Mesh.hpp"
#include "aether/core/ScalarField.hpp"
#include "aether/core/Tensor3x3.hpp"
#include "aether/core/Vector3.hpp"
#include "aether/core/VectorField.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::core;

PYBIND11_MODULE(aether_core_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's C++ math core";

    py::class_<Vector3>(m, "Vector3")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &Vector3::x)
        .def_readwrite("y", &Vector3::y)
        .def_readwrite("z", &Vector3::z)
        .def("dot", &Vector3::dot)
        .def("cross", &Vector3::cross)
        .def("norm", &Vector3::norm)
        .def("normalized", &Vector3::normalized)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(py::self == py::self)
        .def("__repr__", [](const Vector3& v) {
            return "Vector3(" + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " +
                   std::to_string(v.z) + ")";
        });

    py::class_<Tensor3x3>(m, "Tensor3x3")
        .def(py::init<>())
        .def_static("identity", &Tensor3x3::identity)
        .def_static("zero", &Tensor3x3::zero)
        .def("__call__", [](const Tensor3x3& t, int r, int c) { return t(r, c); })
        .def("trace", &Tensor3x3::trace)
        .def("transposed", &Tensor3x3::transposed)
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(py::self * Vector3());

    py::class_<Mesh>(m, "Mesh")
        .def(py::init<>())
        .def("add_vertex", &Mesh::addVertex)
        .def("add_cell", &Mesh::addCell)
        .def("vertex_count", &Mesh::vertexCount)
        .def("cell_count", &Mesh::cellCount)
        .def("vertex", &Mesh::vertex)
        .def("cell", &Mesh::cell)
        .def("cell_centroid", &Mesh::cellCentroid);

    py::class_<ScalarField>(m, "ScalarField")
        .def(py::init<const Mesh&, double>(), py::arg("mesh"), py::arg("initial_value") = 0.0)
        .def("size", &ScalarField::size)
        .def("fill", &ScalarField::fill)
        .def("__getitem__", [](const ScalarField& f, std::size_t i) { return f[i]; })
        .def("__setitem__", [](ScalarField& f, std::size_t i, double v) { f[i] = v; })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double());

    py::class_<VectorField>(m, "VectorField")
        .def(py::init<const Mesh&, const Vector3&>(), py::arg("mesh"), py::arg("initial_value") = Vector3{})
        .def("size", &VectorField::size)
        .def("fill", &VectorField::fill)
        .def("__getitem__", [](const VectorField& f, std::size_t i) { return f[i]; })
        .def("__setitem__", [](VectorField& f, std::size_t i, const Vector3& v) { f[i] = v; })
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double());
}
