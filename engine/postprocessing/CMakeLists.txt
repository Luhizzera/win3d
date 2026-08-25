add_library(aether_postprocessing STATIC
    src/Streamline2D.cpp
    src/MarchingSquares2D.cpp
    src/MarchingCubes3D.cpp
    src/VtkWriter.cpp
)

target_include_directories(aether_postprocessing PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# aether_mesh joins aether_core here because VtkWriter exports a
# TetrahedralMesh: post-processing operates on meshes and fields, so
# depending on the mesh layer is the direction the dependency already
# runs, not a new coupling.
target_link_libraries(aether_postprocessing PUBLIC aether_core aether_mesh)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
