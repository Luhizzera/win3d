if(WIN32)
    add_executable(aether_field_viewer main.cpp)
    target_link_libraries(aether_field_viewer PRIVATE aether_mesh aether_solver aether_app_common)
endif()
