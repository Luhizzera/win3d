#pragma once

#include "aether/plugin/PluginAbi.h"

#include <string>
#include <vector>

namespace aether::plugin {

struct DiagnosticInfo {
    std::string name;
    std::string description;
    std::string pluginName; // which loaded plugin provided it
};

// Module 14 (plugin system), first pass: loads shared libraries at run time
// and exposes the diagnostics they register, so a third party can extend
// the engine with a new field diagnostic without rebuilding it.
//
// The C++ convenience lives entirely on this side of the boundary; the
// boundary itself is pure C (see PluginAbi.h for why that is not
// negotiable). Everything a plugin hands over is copied into owned C++
// types here, so callers never hold a pointer into a plugin's memory.
//
// Loaded libraries stay loaded for this object's lifetime and are released
// in the destructor. **Function pointers obtained from a plugin are only
// valid while that plugin is loaded**, which is why this class owns both
// the handles and the diagnostics rather than handing raw pointers out.
// Non-copyable for the same reason: two hosts unloading the same library
// would be a double-free of the OS handle.
class PluginHost {
public:
    PluginHost() = default;
    ~PluginHost();

    PluginHost(const PluginHost&) = delete;
    PluginHost& operator=(const PluginHost&) = delete;

    // Loads the shared library at `path` and registers its diagnostics.
    // Throws std::runtime_error if the library cannot be opened, does not
    // export AETHER_PLUGIN_REGISTER_SYMBOL, declines to register (returns
    // NULL), or reports an ABI version this build does not support -- the
    // last of which is a refusal to proceed, not a warning, because
    // reading a mismatched struct layout is undefined behavior rather than
    // a recoverable error.
    void load(const std::string& path);

    // Every diagnostic registered by every currently-loaded plugin.
    std::vector<DiagnosticInfo> diagnostics() const;

    bool hasDiagnostic(const std::string& name) const;

    // Runs a plugin-provided diagnostic over `field`. Throws
    // std::invalid_argument if no diagnostic by that name is registered.
    double compute(const std::string& name, const std::vector<double>& field) const;

    std::size_t pluginCount() const { return libraries_.size(); }

private:
    struct LoadedDiagnostic {
        DiagnosticInfo info;
        AetherDiagnosticFn compute;
    };

    // void* rather than HMODULE so this header stays free of <windows.h>
    // (which would leak min/max macros and a few thousand other symbols
    // into every translation unit that includes this). The platform type
    // is cast back inside PluginHost.cpp, the only place that needs it.
    std::vector<void*> libraries_;
    std::vector<LoadedDiagnostic> diagnostics_;
};

} // namespace aether::plugin
