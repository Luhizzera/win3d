add_library(aether_plugin STATIC
    src/PluginHost.cpp
)

target_include_directories(aether_plugin PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# dlopen/dlsym live in libdl on Linux (they are in libc on macOS/BSD, and
# the Windows equivalents need no extra library at all).
if(UNIX AND NOT APPLE)
    target_link_libraries(aether_plugin PRIVATE dl)
endif()

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
