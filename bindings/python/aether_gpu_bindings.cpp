#include "aether/gpu/PoissonOperatorCuda.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace aether::gpu;

// Module 10: GPU. This binding target only exists in the build at all when
// CMake found a CUDA compiler (see bindings/python/CMakeLists.txt) --
// mirrors how every optional piece of this project degrades (pybind11
// itself, at the top level), rather than breaking the rest of the Python
// package on a machine without the CUDA Toolkit.
PYBIND11_MODULE(aether_gpu_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's GPU layer (Module 10, CUDA)";

    py::class_<PoissonOperatorCuda>(m, "PoissonOperatorCuda")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, double, double>(), py::arg("nx"),
             py::arg("ny"), py::arg("nz"), py::arg("dx"), py::arg("dy"), py::arg("dz"))
        .def("apply", &PoissonOperatorCuda::apply, py::arg("x"))
        .def("available", &PoissonOperatorCuda::available);
}
