/* A deliberately incompatible Aether plugin: it is well-formed in every
 * respect except that it reports an ABI version this engine does not
 * support.
 *
 * It exists solely so plugin_tests.cpp can prove the host's version check
 * actually refuses such a plugin -- rather than the test merely asserting
 * that the check exists in the source, which would pass even if the guard
 * were subtly wrong (e.g. comparing the wrong field, or warning instead of
 * refusing). Loading a real mismatched plugin is the only way to
 * demonstrate the refusal path end to end.
 *
 * Its diagnostic deliberately returns an obviously-wrong sentinel: if the
 * host ever *did* load this plugin despite the version mismatch, a test
 * would see -12345.0 rather than a plausible number.
 */

#include "aether/plugin/PluginAbi.h"

#ifdef _WIN32
#define AETHER_PLUGIN_EXPORT __declspec(dllexport)
#else
#define AETHER_PLUGIN_EXPORT __attribute__((visibility("default")))
#endif

static double computeNothing(const double* field, size_t count) {
    (void)field;
    (void)count;
    return -12345.0;
}

static const AetherDiagnostic kDiagnostics[] = {
    {"should_never_load", "If this is ever reachable, the ABI version guard failed", computeNothing},
};

static const AetherPluginInfo kPluginInfo = {
    AETHER_PLUGIN_ABI_VERSION + 998u, /* a version this engine cannot possibly support */
    "bad_version_plugin",
    (uint32_t)(sizeof(kDiagnostics) / sizeof(kDiagnostics[0])),
    kDiagnostics,
};

AETHER_PLUGIN_EXPORT const AetherPluginInfo* aether_plugin_register(void) { return &kPluginInfo; }
