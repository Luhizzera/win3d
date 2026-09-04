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
# OFF only for the wheel build (see pyproject.toml's cmake.define): the
# desktop viewer and its Win32/OpenGL dependency chain are irrelevant to
# `pip install`ing the Python package and cost real build time (and, on a
# platform without a desktop OpenGL stack, could fail outright) for no
# benefit to that build. ON by default so every existing local/CI build
# keeps building the app exactly as before.
option(AETHER_BUILD_APPS "Build the example plugins and desktop viewer app" ON)
# Module 10: GPU. CUDA was chosen deliberately over the dependency-free
# alternatives (GL 4.3 compute shaders, OpenCL) -- see engine/gpu's own
# README/header comments for the tradeoff this reverses (a CUDA Toolkit
# install is now a real external SDK dependency of the C++ build, and the
# engine is tied to NVIDIA hardware for this module specifically). Detected
# rather than required: a machine without the CUDA Toolkit still builds
# every other module untouched, the same optional-dependency pattern
# AETHER_BUILD_PYTHON_BINDINGS already uses for pybind11 below.
option(AETHER_ENABLE_CUDA "Build GPU kernels via CUDA (Module 10)" ON)
# Fase 5: curved-surface STEP/IGES via a real CAD kernel. Detected rather
# than required, the same optional-dependency shape as AETHER_ENABLE_CUDA
# above: a machine without the OpenCASCADE vcpkg toolchain configured still
# builds everything else untouched, including the existing planar-faceted
# STEP loader (engine/geometry/src/StepIO.cpp), which needs no external
# dependency at all. See vcpkg.json (repo root) -- configure with
# -DCMAKE_TOOLCHAIN_FILE=<vcpkg-clone>/scripts/buildsystems/vcpkg.cmake to
# enable it; that toolchain file installs the manifest's own dependencies
# (currently just "opencascade") on first configure.
option(AETHER_ENABLE_OPENCASCADE "Import curved-surface STEP/IGES via OpenCASCADE (Fase 5)" ON)

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

if(AETHER_ENABLE_OPENCASCADE)
    find_package(OpenCASCADE CONFIG QUIET)
    if(NOT OpenCASCADE_FOUND)
        message(STATUS "OpenCASCADE not found - skipping curved-surface STEP/IGES import (Fase 5). "
                        "The planar-faceted STEP loader still builds and works without it. "
                        "Configure with -DCMAKE_TOOLCHAIN_FILE=<vcpkg-clone>/scripts/buildsystems/vcpkg.cmake "
                        "to enable it.")
    endif()
endif()

add_subdirectory(engine/testing)
# engine/gpu added here, before engine/solver: Module 10's own opt-in GPU
# path inside StaggeredCavityBase3D (engine/solver) needs the aether_gpu
# target to already exist by the time engine/solver/CMakeLists.txt runs
# its own `if(TARGET aether_gpu)` check. engine/gpu itself links nothing
# else in this project (confirmed by reading its own CMakeLists.txt), so
# moving it this early changes nothing about what it depends on.
if(AETHER_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
    add_subdirectory(engine/gpu)
endif()
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
# Not gated behind AETHER_BUILD_APPS below: engine/plugin/tests hard-depends
# on these two targets existing whenever AETHER_BUILD_TESTS is on (it passes
# their $<TARGET_FILE:...> in as a compile definition), and they are tiny,
# portable C libraries with no Win32/OpenGL dependency of their own -- the
# thing AETHER_BUILD_APPS actually exists to skip.
add_subdirectory(plugins)
if(AETHER_BUILD_APPS)
    add_subdirectory(apps/common)
    # apps/unified_viewer supersedes the four standalone viewers this
    # project started with (viewer, field_viewer, cavity_viewer,
    # turbulence_viewer): each became a mode of the single executable, so
    # the originals were removed rather than left as a second, drifting
    # copy of the same code.
    add_subdirectory(apps/unified_viewer)
endif()

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
                # Forward slashes, not the backslashes `--cmakedir` prints on
                # Windows: the vcpkg toolchain file (in play once
                # AETHER_ENABLE_OPENCASCADE finds a vcpkg-provided
                # OpenCASCADE) wraps every find_package call through its own
                # macro, which re-parses string arguments and treats a raw
                # backslash as an escape character -- "\J" in
                # "C:\JS\Win3d\..." is not a valid escape, which is a hard
                # configure error, not a warning. file(TO_CMAKE_PATH ...) is
                # the standard normalization for exactly this, and is a
                # no-op for a build without the vcpkg toolchain.
                file(TO_CMAKE_PATH "${AETHER_PYBIND11_DIR}" AETHER_PYBIND11_DIR)
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
