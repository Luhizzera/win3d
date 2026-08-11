add_library(aether_core STATIC
    src/Tensor3x3.cpp
    src/Mesh.cpp
    src/ScalarField.cpp
    src/VectorField.cpp
    src/BigInt.cpp
)

target_include_directories(aether_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
