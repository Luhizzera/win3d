add_executable(aether_plugin_tests plugin_tests.cpp)
target_link_libraries(aether_plugin_tests PRIVATE aether_plugin aether_testing)

# The tests load real shared libraries by path at run time, so they need to
# know where CMake actually put them (which varies by generator: multi-config
# generators like Visual Studio add a per-config subdirectory). A generator
# expression resolves that at build time rather than the test guessing.
target_compile_definitions(aether_plugin_tests PRIVATE
    AETHER_EXAMPLE_PLUGIN_PATH="$<TARGET_FILE:aether_example_plugin>"
    AETHER_BAD_VERSION_PLUGIN_PATH="$<TARGET_FILE:aether_bad_version_plugin>"
)

# Nothing links these, so without an explicit dependency they might not be
# built before the test that loads them runs.
add_dependencies(aether_plugin_tests aether_example_plugin aether_bad_version_plugin)

add_test(NAME aether_plugin_tests COMMAND aether_plugin_tests)
