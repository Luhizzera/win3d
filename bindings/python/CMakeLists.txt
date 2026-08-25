pybind11_add_module(aether_core_py aether_core_bindings.cpp)
target_link_libraries(aether_core_py PRIVATE aether_core)

pybind11_add_module(aether_geometry_py aether_geometry_bindings.cpp)
target_link_libraries(aether_geometry_py PRIVATE aether_geometry)

pybind11_add_module(aether_mesh_py aether_mesh_bindings.cpp)
target_link_libraries(aether_mesh_py PRIVATE aether_mesh)

pybind11_add_module(aether_solver_py aether_solver_bindings.cpp)
target_link_libraries(aether_solver_py PRIVATE aether_solver)

pybind11_add_module(aether_postprocessing_py aether_postprocessing_bindings.cpp)
target_link_libraries(aether_postprocessing_py PRIVATE aether_postprocessing aether_mesh)

pybind11_add_module(aether_persistence_py aether_persistence_bindings.cpp)
target_link_libraries(aether_persistence_py PRIVATE aether_persistence aether_mesh)

pybind11_add_module(aether_optimization_py aether_optimization_bindings.cpp)
target_link_libraries(aether_optimization_py PRIVATE aether_optimization)

pybind11_add_module(aether_analysis_py aether_analysis_bindings.cpp)
target_link_libraries(aether_analysis_py PRIVATE aether_analysis)

pybind11_add_module(aether_ml_py aether_ml_bindings.cpp)
target_link_libraries(aether_ml_py PRIVATE aether_ml)

pybind11_add_module(aether_plugin_py aether_plugin_bindings.cpp)
target_link_libraries(aether_plugin_py PRIVATE aether_plugin)

# aether_gpu only exists as a target when CMake found a CUDA compiler (see
# the root CMakeLists.txt's own check_language(CUDA) guard) -- so this
# binding is conditional the same way, rather than making the CUDA Toolkit
# a hard requirement for the rest of the Python package.
if(TARGET aether_gpu)
    pybind11_add_module(aether_gpu_py aether_gpu_bindings.cpp)
    target_link_libraries(aether_gpu_py PRIVATE aether_gpu)
endif()

# The binding modules were, until now, only ever *compiled* -- CI installs
# pybind11 precisely so they get built, but nothing imported them. That
# catches a signature error and nothing else: a binding can compile and still
# return a dangling pointer, drop a default argument, or mis-order
# constructor parameters. This registers a Python-level suite with ctest so
# the surface is exercised the same way every C++ layer is.
#
# The two arguments are the directory the extension modules were written to
# (a generator expression, since a multi-config generator only knows it at
# build time) and the source tree's Python package, whose own import list is
# part of what is being tested.
if(AETHER_BUILD_TESTS)
    if(DEFINED Python_EXECUTABLE)
        set(aether_python "${Python_EXECUTABLE}")
    elseif(DEFINED PYTHON_EXECUTABLE)
        set(aether_python "${PYTHON_EXECUTABLE}")
    endif()

    if(aether_python)
        add_test(NAME aether_python_bindings_tests
                 COMMAND "${aether_python}"
                         "${CMAKE_SOURCE_DIR}/python/tests/test_unstructured_bindings.py"
                         "$<TARGET_FILE_DIR:aether_solver_py>"
                         "${CMAKE_SOURCE_DIR}/python")
        # Registered separately rather than folded into the suite above:
        # this one exercises the geometry -> volume mesh -> solver pipeline,
        # the only path in the engine that starts from an *imported* surface
        # instead of a synthetic lattice, so a regression there would show
        # up nowhere else. Same argument convention (extension dir, package
        # dir) as the bindings suite.
        add_test(NAME aether_pipeline_tests
                 COMMAND "${aether_python}"
                         "${CMAKE_SOURCE_DIR}/python/tests/test_pipeline.py"
                         "$<TARGET_FILE_DIR:aether_solver_py>"
                         "${CMAKE_SOURCE_DIR}/python")
    else()
        # Printed rather than silently skipped: a suite that quietly stops
        # registering is the exact failure the CI's --no-tests=error flag
        # exists to prevent, and hiding it here would defeat that.
        message(WARNING "Python interpreter not known to CMake - skipping aether_python_bindings_tests")
    endif()
endif()
