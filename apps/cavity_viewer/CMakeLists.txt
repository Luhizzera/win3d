if(WIN32)
    add_executable(aether_cavity_viewer main.cpp)
    target_link_libraries(aether_cavity_viewer PRIVATE aether_solver aether_postprocessing aether_app_common)
endif()
