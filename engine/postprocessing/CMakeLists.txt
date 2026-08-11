add_library(aether_postprocessing STATIC
    src/Streamline2D.cpp
    src/MarchingSquares2D.cpp
    src/MarchingCubes3D.cpp
)

target_include_directories(aether_postprocessing PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(aether_postprocessing PUBLIC aether_core)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
