#ifndef AETHER_PLUGIN_ABI_H
#define AETHER_PLUGIN_ABI_H

/* Module 14 (plugin system): the contract between the Aether engine and a
 * plugin loaded from a shared library at run time.
 *
 * **This header is deliberately pure C, not C++**, and that is the single
 * most important design decision in this module. A plugin is compiled
 * separately from the engine -- potentially by a different compiler, a
 * different version of the same compiler, or the same compiler with
 * different flags (debug vs release CRT, different /std level, different
 * exception model). C++ types have no stable ABI across any of those:
 * passing a std::vector or std::string across a shared-library boundary
 * whose two sides disagree on the layout is undefined behavior that
 * typically manifests as heap corruption rather than a clean error. Only
 * plain C -- POD structs, primitive types, function pointers -- has a
 * stable, documented ABI here. So the boundary is C, and the C++
 * convenience layer (PluginHost) lives entirely on the engine's side of it.
 *
 * Corollaries of that same reasoning, each enforced by convention rather
 * than by the compiler:
 *   - No memory crosses the boundary in a way that requires the other side
 *     to free it. Field data is borrowed (const pointer + count, owned by
 *     the caller for the duration of the call); strings returned by a
 *     plugin must have static storage duration (a string literal, or a
 *     buffer that outlives the plugin's unload).
 *   - No C++ exception may propagate out of a plugin function. Unwinding
 *     across a shared-library boundary between mismatched runtimes is not
 *     defined. A plugin that can fail should report it in its return value.
 *
 * Extension point in this first pass: **named scalar diagnostics over a
 * field** -- a function taking a flat array of doubles and returning one
 * double. Chosen because it is the smallest extension point that is
 * genuinely useful (it is exactly the shape of engine/analysis'
 * FlowDiagnostics functions, so a plugin can add a new diagnostic the
 * engine never shipped) while keeping the ABI trivially safe: nothing is
 * allocated, nothing is freed, nothing is mutated across the boundary.
 */

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Bumped whenever anything below changes shape. The host refuses to load a
 * plugin reporting a different value -- see PluginHost::load(). Without
 * this check, a plugin built against an older layout of AetherPluginInfo
 * would be read using the new layout, i.e. undefined behavior at the very
 * first field access. */
#define AETHER_PLUGIN_ABI_VERSION 1u

/* Computes one scalar from a field of `count` doubles. `field` is borrowed
 * for the duration of the call only and must not be modified or retained.
 * Must not throw; must not return through any path other than a plain
 * return. */
typedef double (*AetherDiagnosticFn)(const double* field, size_t count);

typedef struct AetherDiagnostic {
    const char* name;        /* static storage duration; identifies this diagnostic to callers */
    const char* description; /* static storage duration; human-readable, may be shown in a UI */
    AetherDiagnosticFn compute;
} AetherDiagnostic;

typedef struct AetherPluginInfo {
    uint32_t abiVersion; /* must equal AETHER_PLUGIN_ABI_VERSION */
    const char* pluginName;
    uint32_t diagnosticCount;
    const AetherDiagnostic* diagnostics; /* array of diagnosticCount entries, static storage duration */
} AetherPluginInfo;

/* The one symbol every Aether plugin must export, by exactly this name.
 * Returns a pointer with static storage duration (the engine does not free
 * it) describing everything the plugin provides. Returning NULL means the
 * plugin declines to load. */
typedef const AetherPluginInfo* (*AetherPluginRegisterFn)(void);

#define AETHER_PLUGIN_REGISTER_SYMBOL "aether_plugin_register"

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* AETHER_PLUGIN_ABI_H */
