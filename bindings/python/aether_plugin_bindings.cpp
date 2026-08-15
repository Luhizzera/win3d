#include "aether/plugin/PluginHost.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::plugin;

PYBIND11_MODULE(aether_plugin_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's plugin system (Module 14)";

    m.attr("ABI_VERSION") = AETHER_PLUGIN_ABI_VERSION;

    py::class_<DiagnosticInfo>(m, "DiagnosticInfo")
        .def_readonly("name", &DiagnosticInfo::name)
        .def_readonly("description", &DiagnosticInfo::description)
        .def_readonly("plugin_name", &DiagnosticInfo::pluginName);

    // PluginHost owns the loaded libraries and unloads them on destruction,
    // and is non-copyable for exactly that reason -- so it is bound without
    // any copy support. Python's own reference counting then governs when
    // the libraries are released, which is the behavior a caller expects.
    py::class_<PluginHost>(m, "PluginHost")
        .def(py::init<>())
        .def("load", &PluginHost::load, py::arg("path"))
        .def("diagnostics", &PluginHost::diagnostics)
        .def("has_diagnostic", &PluginHost::hasDiagnostic, py::arg("name"))
        .def("compute", &PluginHost::compute, py::arg("name"), py::arg("field"))
        .def("plugin_count", &PluginHost::pluginCount);
}
