/* An example Aether plugin, deliberately written in **plain C** rather than
 * C++ -- proving the point PluginAbi.h makes: the plugin boundary is a C
 * ABI, so a plugin needs no C++ runtime compatibility with the engine at
 * all, and can be written in any language able to export a C symbol.
 *
 * It lives outside engine/ on purpose: this is what a third party's plugin
 * would look like, and it is built as a separate shared library that the
 * engine knows nothing about until it is loaded by path at run time.
 *
 * The two diagnostics below were chosen because each has an exact,
 * independently-computable closed form, so the test that loads this plugin
 * can verify the value that came back across the library boundary against
 * the same quantity computed in-process -- see plugin_tests.cpp.
 */

#include "aether/plugin/PluginAbi.h"

#include <math.h>

#ifdef _WIN32
#define AETHER_PLUGIN_EXPORT __declspec(dllexport)
#else
#define AETHER_PLUGIN_EXPORT __attribute__((visibility("default")))
#endif

/* Root mean square: sqrt(sum(x_i^2) / n). */
static double computeRms(const double* field, size_t count) {
    double sumOfSquares = 0.0;
    size_t i;
    if (count == 0) {
        return 0.0;
    }
    for (i = 0; i < count; ++i) {
        sumOfSquares += field[i] * field[i];
    }
    return sqrt(sumOfSquares / (double)count);
}

/* Peak-to-peak amplitude: max - min. */
static double computePeakToPeak(const double* field, size_t count) {
    double minValue;
    double maxValue;
    size_t i;
    if (count == 0) {
        return 0.0;
    }
    minValue = field[0];
    maxValue = field[0];
    for (i = 1; i < count; ++i) {
        if (field[i] < minValue) {
            minValue = field[i];
        }
        if (field[i] > maxValue) {
            maxValue = field[i];
        }
    }
    return maxValue - minValue;
}

/* Static storage duration, as the ABI requires: these outlive every call
 * and are never freed by the engine. */
static const AetherDiagnostic kDiagnostics[] = {
    {"rms", "Root mean square of the field", computeRms},
    {"peak_to_peak", "Maximum minus minimum value of the field", computePeakToPeak},
};

static const AetherPluginInfo kPluginInfo = {
    AETHER_PLUGIN_ABI_VERSION,
    "example_diagnostics",
    (uint32_t)(sizeof(kDiagnostics) / sizeof(kDiagnostics[0])),
    kDiagnostics,
};

AETHER_PLUGIN_EXPORT const AetherPluginInfo* aether_plugin_register(void) { return &kPluginInfo; }
