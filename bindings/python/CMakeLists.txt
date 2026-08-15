pybind11_add_module(aether_core_py aether_core_bindings.cpp)
target_link_libraries(aether_core_py PRIVATE aether_core)

pybind11_add_module(aether_geometry_py aether_geometry_bindings.cpp)
target_link_libraries(aether_geometry_py PRIVATE aether_geometry)

pybind11_add_module(aether_mesh_py aether_mesh_bindings.cpp)
target_link_libraries(aether_mesh_py PRIVATE aether_mesh)

pybind11_add_module(aether_solver_py aether_solver_bindings.cpp)
target_link_libraries(aether_solver_py PRIVATE aether_solver)

pybind11_add_module(aether_postprocessing_py aether_postprocessing_bindings.cpp)
target_link_libraries(aether_postprocessing_py PRIVATE aether_postprocessing)

pybind11_add_module(aether_persistence_py aether_persistence_bindings.cpp)
target_link_libraries(aether_persistence_py PRIVATE aether_persistence)

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
