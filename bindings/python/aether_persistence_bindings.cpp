#include "aether/persistence/FieldArchive.hpp"
#include "aether/persistence/GridArchive.hpp"
#include "aether/persistence/ProjectHistory.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::persistence;

PYBIND11_MODULE(aether_persistence_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's persistence layer (Module 11)";

    py::class_<FieldArchive>(m, "FieldArchive")
        .def(py::init<>())
        .def("set_metadata", &FieldArchive::setMetadata, py::arg("key"), py::arg("value"))
        .def("metadata", &FieldArchive::metadata, py::arg("key"))
        .def("has_metadata", &FieldArchive::hasMetadata, py::arg("key"))
        .def("set_field", &FieldArchive::setField, py::arg("name"), py::arg("data"))
        .def("field", &FieldArchive::field, py::arg("name"))
        .def("has_field", &FieldArchive::hasField, py::arg("name"))
        .def("save", &FieldArchive::save, py::arg("path"))
        .def_static("load", &FieldArchive::load, py::arg("path"));

    py::class_<HistoryEntry>(m, "HistoryEntry")
        .def_readonly("label", &HistoryEntry::label)
        .def_readonly("filename", &HistoryEntry::filename)
        .def_readonly("timestamp", &HistoryEntry::timestamp);

    py::class_<ProjectHistory>(m, "ProjectHistory")
        .def(py::init<std::string>(), py::arg("directory"))
        .def("record", &ProjectHistory::record, py::arg("label"), py::arg("archive"))
        .def("entries", &ProjectHistory::entries)
        .def("load", &ProjectHistory::load, py::arg("entry"));

    m.def("save_grid", &saveGrid, py::arg("archive"), py::arg("grid"));
    m.def("load_grid", &loadGrid, py::arg("archive"));
}
