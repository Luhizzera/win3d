add_library(aether_geometry STATIC
    src/TriangleMesh.cpp
    src/StlIO.cpp
    src/ObjIO.cpp
    src/StepIO.cpp
)

target_include_directories(aether_geometry PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# aether_mesh (not just aether_core) because StepIO reuses
# PolygonTriangulation2D to triangulate each planar STEP face -- no
# circular dependency, aether_mesh only depends on aether_core itself.
target_link_libraries(aether_geometry PUBLIC aether_core aether_mesh)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
