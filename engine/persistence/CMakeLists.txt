add_library(aether_persistence STATIC
    src/FieldArchive.cpp
    src/ProjectHistory.cpp
    src/GridArchive.cpp
)

target_include_directories(aether_persistence PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(aether_persistence PUBLIC aether_mesh)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
