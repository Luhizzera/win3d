# Example plugins, built as real shared libraries loaded by path at run
# time -- the engine links against neither of them. They sit outside
# engine/ deliberately: this is the shape a third party's plugin would
# take.
#
# Both are C, not C++, which is the point (see PluginAbi.h): a plugin needs
# only a C compiler and this project's one C header, no C++ runtime
# compatibility with the engine whatsoever.

function(aether_add_plugin name source)
    add_library(${name} SHARED ${source})
    target_include_directories(${name} PRIVATE ${CMAKE_SOURCE_DIR}/engine/plugin/include)
    set_target_properties(${name} PROPERTIES C_STANDARD 99 C_STANDARD_REQUIRED ON)
endfunction()

aether_add_plugin(aether_example_plugin example_diagnostics/example_diagnostics.c)

# Deliberately reports an unsupported ABI version; exists so the test suite
# can prove PluginHost actually refuses such a plugin. See the source file.
aether_add_plugin(aether_bad_version_plugin bad_version_plugin/bad_version_plugin.c)
