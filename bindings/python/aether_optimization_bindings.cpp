#include "aether/optimization/NelderMead.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::optimization;

PYBIND11_MODULE(aether_optimization_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's optimization layer (Module 12)";

    py::class_<OptimizationResult>(m, "OptimizationResult")
        .def_readonly("parameters", &OptimizationResult::parameters)
        .def_readonly("value", &OptimizationResult::value)
        .def_readonly("iterations", &OptimizationResult::iterations)
        .def_readonly("converged", &OptimizationResult::converged);

    py::class_<NelderMead>(m, "NelderMead")
        .def(py::init<std::size_t, double>(), py::arg("max_iterations") = 1000, py::arg("tolerance") = 1e-10)
        .def("minimize", &NelderMead::minimize, py::arg("objective"), py::arg("initial_guess"),
             py::arg("initial_step_size") = 0.1);
}
