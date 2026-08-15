#include "aether/ml/MultiLayerPerceptron.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::ml;

PYBIND11_MODULE(aether_ml_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's surrogate-model layer (Module 12)";

    py::class_<MultiLayerPerceptron>(m, "MultiLayerPerceptron")
        .def(py::init<const std::vector<std::size_t>&, unsigned>(), py::arg("layer_sizes"), py::arg("seed") = 42)
        .def("predict", &MultiLayerPerceptron::predict, py::arg("input"))
        .def("loss", &MultiLayerPerceptron::loss, py::arg("inputs"), py::arg("targets"))
        .def("gradient", &MultiLayerPerceptron::gradient, py::arg("inputs"), py::arg("targets"))
        .def("parameters", &MultiLayerPerceptron::parameters)
        .def("set_parameters", &MultiLayerPerceptron::setParameters, py::arg("params"))
        .def("train_epoch", &MultiLayerPerceptron::trainEpoch, py::arg("inputs"), py::arg("targets"),
             py::arg("learning_rate"));
}
