add_library(aether_ml STATIC
    src/MultiLayerPerceptron.cpp
)

target_include_directories(aether_ml PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
