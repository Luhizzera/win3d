if(WIN32)
    add_executable(aether_turbulence_viewer main.cpp)
    target_link_libraries(aether_turbulence_viewer PRIVATE aether_solver aether_app_common)
endif()
