add_library(aether_optimization STATIC
    src/NelderMead.cpp
)

target_include_directories(aether_optimization PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
