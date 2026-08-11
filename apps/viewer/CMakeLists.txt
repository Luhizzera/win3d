if(WIN32)
    add_executable(aether_viewer main.cpp)
    target_link_libraries(aether_viewer PRIVATE aether_core aether_geometry aether_app_common)
endif()
