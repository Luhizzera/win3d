#include "aether/analysis/FlowDiagnostics.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::analysis;

PYBIND11_MODULE(aether_analysis_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's analytical diagnostics layer (Module 12)";

    py::class_<FieldStatistics>(m, "FieldStatistics")
        .def_readonly("min_value", &FieldStatistics::minValue)
        .def_readonly("max_value", &FieldStatistics::maxValue)
        .def_readonly("mean", &FieldStatistics::mean);

    m.def("compute_statistics", &computeStatistics, py::arg("field"));
    m.def("max_courant_number", &maxCourantNumber, py::arg("u"), py::arg("v"), py::arg("dx"), py::arg("dy"),
          py::arg("dt"));
    m.def("checkerboard_index", &checkerboardIndex, py::arg("field"), py::arg("nx"), py::arg("ny"));
    m.def("summarize_field", &summarizeField, py::arg("name"), py::arg("field"));
}
