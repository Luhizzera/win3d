cmake_minimum_required(VERSION 3.20)
project(Aether LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(AETHER_BUILD_TESTS "Build Aether C++ unit tests" ON)
option(AETHER_BUILD_PYTHON_BINDINGS "Build Python bindings (requires pybind11)" ON)

enable_testing()

add_subdirectory(engine/testing)
add_subdirectory(engine/core)
add_subdirectory(engine/geometry)
add_subdirectory(engine/mesh)
add_subdirectory(engine/solver)
add_subdirectory(engine/postprocessing)
add_subdirectory(apps/common)
# apps/unified_viewer supersedes the four standalone viewers this project
# started with (viewer, field_viewer, cavity_viewer, turbulence_viewer):
# each became a mode of the single executable, so the originals were
# removed rather than left as a second, drifting copy of the same code.
add_subdirectory(apps/unified_viewer)

if(AETHER_BUILD_PYTHON_BINDINGS)
    find_package(pybind11 CONFIG QUIET)
    if(pybind11_FOUND)
        add_subdirectory(bindings/python)
    else()
        message(STATUS "pybind11 not found - skipping Python bindings. "
                        "Run 'pip install pybind11' in your Python env and re-run cmake to enable them.")
    endif()
endif()
