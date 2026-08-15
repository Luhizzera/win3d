add_library(aether_analysis STATIC
    src/FlowDiagnostics.cpp
)

target_include_directories(aether_analysis PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
