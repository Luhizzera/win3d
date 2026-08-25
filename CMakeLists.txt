cmake_minimum_required(VERSION 3.20)
# C is enabled alongside C++ for Module 14: the plugin ABI is pure C (see
# engine/plugin/include/aether/plugin/PluginAbi.h for why that is not
# negotiable), and the example plugins under plugins/ are written in C to
# demonstrate that a plugin needs no C++ runtime compatibility with the
# engine at all. No engine library itself is C.
project(Aether LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(AETHER_BUILD_TESTS "Build Aether C++ unit tests" ON)
option(AETHER_BUILD_PYTHON_BINDINGS "Build Python bindings (requires pybind11)" ON)
# Module 10: GPU. CUDA was chosen deliberately over the dependency-free
# alternatives (GL 4.3 compute shaders, OpenCL) -- see engine/gpu's own
# README/header comments for the tradeoff this reverses (a CUDA Toolkit
# install is now a real external SDK dependency of the C++ build, and the
# engine is tied to NVIDIA hardware for this module specifically). Detected
# rather than required: a machine without the CUDA Toolkit still builds
# every other module untouched, the same optional-dependency pattern
# AETHER_BUILD_PYTHON_BINDINGS already uses for pybind11 below.
option(AETHER_ENABLE_CUDA "Build GPU kernels via CUDA (Module 10)" ON)

enable_testing()

if(AETHER_ENABLE_CUDA)
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        set(CMAKE_CUDA_STANDARD 17)
        set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    else()
        message(STATUS "CUDA compiler not found - skipping engine/gpu (Module 10). "
                        "Install the CUDA Toolkit and re-run cmake to enable it.")
    endif()
endif()

add_subdirectory(engine/testing)
add_subdirectory(engine/core)
add_subdirectory(engine/geometry)
add_subdirectory(engine/mesh)
add_subdirectory(engine/solver)
add_subdirectory(engine/postprocessing)
add_subdirectory(engine/persistence)
add_subdirectory(engine/optimization)
add_subdirectory(engine/analysis)
add_subdirectory(engine/ml)
add_subdirectory(engine/plugin)
# Example plugins (Module 14). Built as standalone shared libraries that
# nothing links against -- the engine discovers them by path at run time.
add_subdirectory(plugins)
if(AETHER_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
    add_subdirectory(engine/gpu)
endif()
add_subdirectory(apps/common)
# apps/unified_viewer supersedes the four standalone viewers this project
# started with (viewer, field_viewer, cavity_viewer, turbulence_viewer):
# each became a mode of the single executable, so the originals were
# removed rather than left as a second, drifting copy of the same code.
add_subdirectory(apps/unified_viewer)

if(AETHER_BUILD_PYTHON_BINDINGS)
    find_package(pybind11 CONFIG QUIET)

    # **Ask a Python interpreter where pybind11 is, before giving up.**
    # find_package alone searches CMake's own prefixes, which do not include
    # a virtualenv's site-packages -- so on a machine where pybind11 is
    # installed *in the venv this repository ships*, the plain lookup fails
    # and the whole Python layer is skipped.
    #
    # That failure was silent and expensive: a clean configure produced a
    # C++-only build, ctest reported 11 suites passing, and the two Python
    # suites simply were not there to fail. `--no-tests=error` does not
    # catch it either, because tests do exist -- just fewer of them. It went
    # unnoticed for as long as it did because an older build directory
    # carried a cached pybind11_DIR from a one-off manual configure, so only
    # a genuinely clean build (or a fresh clone) ever hit it.
    if(NOT pybind11_FOUND)
        find_program(AETHER_PYTHON
                     NAMES python python3
                     HINTS "${CMAKE_SOURCE_DIR}/python/.venv/Scripts"
                           "${CMAKE_SOURCE_DIR}/python/.venv/bin"
                     NO_CACHE)
        if(AETHER_PYTHON)
            execute_process(COMMAND "${AETHER_PYTHON}" -m pybind11 --cmakedir
                            OUTPUT_VARIABLE AETHER_PYBIND11_DIR
                            OUTPUT_STRIP_TRAILING_WHITESPACE
                            ERROR_QUIET)
            if(AETHER_PYBIND11_DIR)
                find_package(pybind11 CONFIG QUIET
                             PATHS "${AETHER_PYBIND11_DIR}" NO_DEFAULT_PATH)
            endif()
        endif()
    endif()

    if(pybind11_FOUND)
        add_subdirectory(bindings/python)
    else()
        # A warning, not a status message: skipping this drops two whole
        # test suites and the entire orchestration layer, which is not a
        # detail to scroll past.
        message(WARNING "pybind11 not found - skipping Python bindings AND the two Python "
                        "test suites that go with them. Run 'pip install pybind11' in the "
                        "Python environment CMake is using, or pass -Dpybind11_DIR=<path>.")
    endif()
endif()
