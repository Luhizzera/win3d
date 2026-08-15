add_library(aether_mesh STATIC
    src/StructuredGrid3D.cpp
    src/DelaunayTriangulation2D.cpp
    src/PolygonTriangulation2D.cpp
    src/DelaunayTetrahedralization3D.cpp
    src/RobustPredicates.cpp
    src/TetrahedralMesh.cpp
)

target_include_directories(aether_mesh PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(aether_mesh PUBLIC aether_core)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
