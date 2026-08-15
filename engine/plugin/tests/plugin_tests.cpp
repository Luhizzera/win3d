#include "aether/plugin/PluginHost.hpp"
#include "aether/testing/Check.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

using aether::plugin::DiagnosticInfo;
using aether::plugin::PluginHost;

namespace {

// A deterministic, non-trivial field: no repeated values, both signs, and
// no symmetry that could make a wrong implementation accidentally agree
// with the right one.
std::vector<double> testField() {
    return {1.5, -2.25, 0.75, 4.0, -0.5, 3.125, -1.0, 2.5};
}

// Independent reimplementations of what the example plugin computes,
// written here from each quantity's own definition rather than shared with
// the plugin's source -- the same "independent reimplementation as
// cross-check" discipline used in aether_gpu_tests and the mesh module.
// If the plugin and the engine ever disagree, that shows up as a real
// mismatch instead of being hidden by shared code.
double referenceRms(const std::vector<double>& field) {
    double sumOfSquares = 0.0;
    for (double value : field) {
        sumOfSquares += value * value;
    }
    return std::sqrt(sumOfSquares / static_cast<double>(field.size()));
}

double referencePeakToPeak(const std::vector<double>& field) {
    const auto [minIt, maxIt] = std::minmax_element(field.begin(), field.end());
    return *maxIt - *minIt;
}

// The core claim of this module: a shared library the engine does not link
// against, built by a different compiler front-end (C, not C++), can be
// loaded at run time and its functions called through the C ABI with
// results identical to computing the same thing in-process.
//
// Equality is expected to be **exact**, not approximate: the plugin and the
// reference below perform the same operations on the same binary64 values
// in the same order, and crossing a library boundary does not change
// arithmetic. A difference here would mean the two implementations really
// are computing different things (or that data was corrupted crossing the
// boundary), not floating-point noise to be tolerated.
void testLoadsPluginAndComputesExactly() {
    PluginHost host;
    host.load(AETHER_EXAMPLE_PLUGIN_PATH);

    AETHER_CHECK(host.pluginCount() == 1);

    const std::vector<DiagnosticInfo> available = host.diagnostics();
    AETHER_CHECK(available.size() == 2);
    for (const DiagnosticInfo& info : available) {
        std::printf("  [aether_plugin_tests] registered: %s (%s) from plugin '%s'\n", info.name.c_str(),
                    info.description.c_str(), info.pluginName.c_str());
        AETHER_CHECK(info.pluginName == "example_diagnostics");
        AETHER_CHECK(!info.description.empty());
    }

    AETHER_CHECK(host.hasDiagnostic("rms"));
    AETHER_CHECK(host.hasDiagnostic("peak_to_peak"));
    AETHER_CHECK(!host.hasDiagnostic("not_a_real_diagnostic"));

    const std::vector<double> field = testField();

    const double pluginRms = host.compute("rms", field);
    const double expectedRms = referenceRms(field);
    std::printf("  [aether_plugin_tests] rms: plugin=%.17g reference=%.17g\n", pluginRms, expectedRms);
    AETHER_CHECK(pluginRms == expectedRms);

    const double pluginRange = host.compute("peak_to_peak", field);
    const double expectedRange = referencePeakToPeak(field);
    std::printf("  [aether_plugin_tests] peak_to_peak: plugin=%.17g reference=%.17g\n", pluginRange,
                expectedRange);
    AETHER_CHECK(pluginRange == expectedRange);

    // Guards against the whole test passing trivially on a degenerate
    // field where every diagnostic would return 0 regardless.
    AETHER_CHECK(pluginRms > 0.0);
    AETHER_CHECK(pluginRange > 0.0);
}

// Proves the ABI version guard by loading a plugin that is well-formed in
// every way *except* its reported version -- the only way to exercise the
// refusal path for real. Asserting the check exists in the source would
// pass even if the guard compared the wrong field or merely warned.
void testRefusesIncompatibleAbiVersion() {
    PluginHost host;
    bool threw = false;
    std::string message;
    try {
        host.load(AETHER_BAD_VERSION_PLUGIN_PATH);
    } catch (const std::runtime_error& exc) {
        threw = true;
        message = exc.what();
    }
    std::printf("  [aether_plugin_tests] bad-version plugin refused with: %s\n", message.c_str());
    AETHER_CHECK(threw);
    AETHER_CHECK(message.find("ABI version") != std::string::npos);

    // The refusal must be complete, not partial: nothing from the rejected
    // plugin may remain registered.
    AETHER_CHECK(host.pluginCount() == 0);
    AETHER_CHECK(host.diagnostics().empty());
    AETHER_CHECK(!host.hasDiagnostic("should_never_load"));
}

void testRejectsMissingAndNonPluginLibraries() {
    PluginHost host;

    bool missingThrew = false;
    std::string missingMessage;
    try {
        host.load("definitely_not_a_real_library_path.dll");
    } catch (const std::runtime_error& exc) {
        missingThrew = true;
        missingMessage = exc.what();
    }
    AETHER_CHECK(missingThrew);
    AETHER_CHECK(missingMessage.find("could not open") != std::string::npos);

    // A genuinely different failure from the one above, and worth
    // separating: a library that opens perfectly well but simply is not an
    // Aether plugin. A well-known system library is the honest way to test
    // that path -- using another made-up filename would just re-test
    // "could not open" while claiming to test something else.
#ifdef _WIN32
    const char* realNonPlugin = "kernel32.dll";
#elif defined(__APPLE__)
    const char* realNonPlugin = "/usr/lib/libSystem.B.dylib";
#else
    const char* realNonPlugin = "libm.so.6";
#endif
    bool notAPluginThrew = false;
    std::string notAPluginMessage;
    try {
        host.load(realNonPlugin);
    } catch (const std::runtime_error& exc) {
        notAPluginThrew = true;
        notAPluginMessage = exc.what();
    }
    std::printf("  [aether_plugin_tests] non-plugin library refused with: %s\n", notAPluginMessage.c_str());
    AETHER_CHECK(notAPluginThrew);
    AETHER_CHECK(notAPluginMessage.find("does not export") != std::string::npos);

    // Neither failure may leave anything behind.
    AETHER_CHECK(host.pluginCount() == 0);
    AETHER_CHECK(host.diagnostics().empty());
}

void testComputeRejectsUnknownDiagnostic() {
    PluginHost host;
    host.load(AETHER_EXAMPLE_PLUGIN_PATH);
    bool threw = false;
    try {
        host.compute("no_such_diagnostic", {1.0, 2.0});
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    AETHER_CHECK(threw);
}

// Two hosts loading the same library, and a host loading two plugins, are
// both ordinary situations a plugin system must survive -- worth checking
// explicitly because the destructor unloads libraries, and getting the
// ownership wrong here would show up as a crash at scope exit rather than
// a failed assertion.
void testMultiplePluginsAndHostsCoexist() {
    {
        PluginHost first;
        PluginHost second;
        first.load(AETHER_EXAMPLE_PLUGIN_PATH);
        second.load(AETHER_EXAMPLE_PLUGIN_PATH);
        AETHER_CHECK(first.compute("rms", testField()) == second.compute("rms", testField()));
    }

    PluginHost host;
    host.load(AETHER_EXAMPLE_PLUGIN_PATH);
    host.load(AETHER_EXAMPLE_PLUGIN_PATH); // same library twice
    AETHER_CHECK(host.pluginCount() == 2);
    AETHER_CHECK(host.diagnostics().size() == 4);
}

} // namespace

int main() {
    testLoadsPluginAndComputesExactly();
    testRefusesIncompatibleAbiVersion();
    testRejectsMissingAndNonPluginLibraries();
    testComputeRejectsUnknownDiagnostic();
    testMultiplePluginsAndHostsCoexist();
    std::printf("aether_plugin_tests: OK\n");
    return 0;
}
