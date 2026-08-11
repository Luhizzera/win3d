add_library(aether_geometry STATIC
    src/TriangleMesh.cpp
    src/StlIO.cpp
    src/ObjIO.cpp
)

target_include_directories(aether_geometry PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(aether_geometry PUBLIC aether_core)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
