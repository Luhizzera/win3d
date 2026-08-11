if(WIN32)
    add_executable(aether_unified_viewer main.cpp)
    target_link_libraries(aether_unified_viewer PRIVATE aether_core aether_geometry aether_mesh aether_solver
                                                          aether_postprocessing aether_app_common)
endif()
