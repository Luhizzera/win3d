#include "aether/plugin/PluginHost.hpp"

#include <stdexcept>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#else
#include <dlfcn.h>
#endif

namespace aether::plugin {

namespace {

// The two platform primitives this module needs, isolated here so the rest
// of the file reads the same on either OS. This project's *apps* are Win32
// only, but every engine/ library so far has been portable C++, and the
// POSIX path costs three lines -- so keeping it portable is cheaper than
// making engine/plugin the first Windows-locked engine library.
void* openLibrary(const std::string& path) {
#ifdef _WIN32
    return static_cast<void*>(::LoadLibraryA(path.c_str()));
#else
    return ::dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
#endif
}

void* findSymbol(void* handle, const char* name) {
#ifdef _WIN32
    return reinterpret_cast<void*>(::GetProcAddress(static_cast<HMODULE>(handle), name));
#else
    return ::dlsym(handle, name);
#endif
}

void closeLibrary(void* handle) {
#ifdef _WIN32
    ::FreeLibrary(static_cast<HMODULE>(handle));
#else
    ::dlclose(handle);
#endif
}

} // namespace

PluginHost::~PluginHost() {
    // Diagnostics point into plugin memory, so they must be dropped before
    // the libraries backing them are unloaded, not after.
    diagnostics_.clear();
    for (void* handle : libraries_) {
        closeLibrary(handle);
    }
}

void PluginHost::load(const std::string& path) {
    void* handle = openLibrary(path);
    if (handle == nullptr) {
        throw std::runtime_error("PluginHost::load: could not open shared library '" + path + "'");
    }

    void* symbol = findSymbol(handle, AETHER_PLUGIN_REGISTER_SYMBOL);
    if (symbol == nullptr) {
        closeLibrary(handle);
        throw std::runtime_error("PluginHost::load: '" + path + "' does not export " +
                                  AETHER_PLUGIN_REGISTER_SYMBOL + " -- not an Aether plugin");
    }

    const auto registerFn = reinterpret_cast<AetherPluginRegisterFn>(symbol);
    const AetherPluginInfo* info = registerFn();
    if (info == nullptr) {
        closeLibrary(handle);
        throw std::runtime_error("PluginHost::load: '" + path + "' declined to register (returned null)");
    }

    // Checked before *any* other field of `info` is read: if the version
    // disagrees, the struct's layout may disagree too, so every subsequent
    // access would be undefined behavior. This is why the check is a hard
    // refusal rather than a warning-and-continue.
    if (info->abiVersion != AETHER_PLUGIN_ABI_VERSION) {
        const std::string reported = std::to_string(info->abiVersion);
        closeLibrary(handle);
        throw std::runtime_error("PluginHost::load: '" + path + "' reports plugin ABI version " + reported +
                                  ", but this build supports version " +
                                  std::to_string(AETHER_PLUGIN_ABI_VERSION) + " -- refusing to load");
    }

    const std::string pluginName = info->pluginName != nullptr ? info->pluginName : "(unnamed)";
    for (std::uint32_t i = 0; i < info->diagnosticCount; ++i) {
        const AetherDiagnostic& entry = info->diagnostics[i];
        if (entry.name == nullptr || entry.compute == nullptr) {
            continue; // a malformed entry is skipped rather than taking down the whole plugin
        }
        LoadedDiagnostic loaded;
        loaded.info.name = entry.name;
        loaded.info.description = entry.description != nullptr ? entry.description : "";
        loaded.info.pluginName = pluginName;
        loaded.compute = entry.compute;
        diagnostics_.push_back(std::move(loaded));
    }

    libraries_.push_back(handle);
}

std::vector<DiagnosticInfo> PluginHost::diagnostics() const {
    std::vector<DiagnosticInfo> result;
    result.reserve(diagnostics_.size());
    for (const LoadedDiagnostic& entry : diagnostics_) {
        result.push_back(entry.info);
    }
    return result;
}

bool PluginHost::hasDiagnostic(const std::string& name) const {
    for (const LoadedDiagnostic& entry : diagnostics_) {
        if (entry.info.name == name) {
            return true;
        }
    }
    return false;
}

double PluginHost::compute(const std::string& name, const std::vector<double>& field) const {
    for (const LoadedDiagnostic& entry : diagnostics_) {
        if (entry.info.name == name) {
            return entry.compute(field.data(), field.size());
        }
    }
    throw std::invalid_argument("PluginHost::compute: no loaded plugin provides a diagnostic named '" + name + "'");
}

} // namespace aether::plugin
