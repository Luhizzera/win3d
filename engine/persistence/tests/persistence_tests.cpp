#include "aether/persistence/FieldArchive.hpp"
#include "aether/persistence/GridArchive.hpp"
#include "aether/persistence/ProjectHistory.hpp"
#include "aether/solver/LidDrivenCavitySolver2D.hpp"
#include "aether/solver/StaggeredLidDrivenCavitySolver3D.hpp"
#include "aether/testing/Check.hpp"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <vector>

using aether::core::Vector3;
using aether::mesh::StructuredGrid3D;
using aether::persistence::FieldArchive;
using aether::persistence::HistoryEntry;
using aether::persistence::ProjectHistory;
using aether::solver::LidDrivenCavitySolver2D;
using aether::solver::StaggeredLidDrivenCavitySolver3D;

namespace {

// Generic round trip: arbitrary metadata + fields, no solver involved.
// Every byte a FieldArchive writes is either raw metadata doubles or raw
// field doubles copied verbatim from a std::vector<double> -- no lossy
// transformation happens anywhere in save()/load(), so equality here is
// expected to be exact, not approximate.
void testGenericRoundTripIsExact() {
    FieldArchive archive;
    archive.setMetadata("nx", 4.0);
    archive.setMetadata("viscosity", 0.01);
    archive.setMetadata("time", 1.25);
    archive.setField("u", {0.1, -0.2, 0.3, 0.0, 5.5});
    archive.setField("p", {1.0, 2.0, 3.0});
    archive.setField("empty", {});

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_persistence_test_generic.aecf";
    archive.save(path.string());
    const FieldArchive loaded = FieldArchive::load(path.string());
    std::filesystem::remove(path);

    AETHER_CHECK(loaded.metadata("nx") == 4.0);
    AETHER_CHECK(loaded.metadata("viscosity") == 0.01);
    AETHER_CHECK(loaded.metadata("time") == 1.25);
    AETHER_CHECK(!loaded.hasMetadata("does_not_exist"));

    const std::vector<double> expectedU = {0.1, -0.2, 0.3, 0.0, 5.5};
    AETHER_CHECK(loaded.field("u") == expectedU);
    const std::vector<double> expectedP = {1.0, 2.0, 3.0};
    AETHER_CHECK(loaded.field("p") == expectedP);
    AETHER_CHECK(loaded.hasField("empty"));
    AETHER_CHECK(loaded.field("empty").empty());
    AETHER_CHECK(!loaded.hasField("does_not_exist"));
}

// Loading a file that isn't an Aether checkpoint must fail loudly (wrong
// magic), not silently return a garbage-filled archive.
void testLoadRejectsBadMagic() {
    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_persistence_test_bad_magic.aecf";
    {
        std::ofstream out(path, std::ios::binary);
        out << "not an archive";
    }

    bool threw = false;
    try {
        FieldArchive::load(path.string());
    } catch (const std::runtime_error&) {
        threw = true;
    }
    std::filesystem::remove(path);
    AETHER_CHECK(threw);
}

// Requesting a metadata key or field name that was never set must fail
// loudly too (std::map::at throws std::out_of_range), not default to 0 /
// an empty vector and hide a caller bug.
void testMissingKeysThrow() {
    FieldArchive archive;
    bool metadataThrew = false;
    try {
        archive.metadata("missing");
    } catch (const std::out_of_range&) {
        metadataThrew = true;
    }
    AETHER_CHECK(metadataThrew);

    bool fieldThrew = false;
    try {
        archive.field("missing");
    } catch (const std::out_of_range&) {
        fieldThrew = true;
    }
    AETHER_CHECK(fieldThrew);
}

// The actual point of Module 11's first pass: a real solver's state,
// extracted field by field into a FieldArchive, survives a save/load round
// trip bit-exactly. This does not yet resume a *live* solver from a
// checkpoint (LidDrivenCavitySolver2D has no setter to inject external
// field data back in) -- that is deliberately left as follow-up work, the
// same way Module 10 first validated PoissonOperatorCuda in isolation
// before any solver used it. What this proves is that the archive format
// itself is a faithful, lossless container for genuine CFD field data, not
// just for hand-picked test vectors.
void testRoundTripsRealSolverState() {
    LidDrivenCavitySolver2D solver(8, 6, 1.0, 0.75, 0.01, 1.0);
    for (int step = 0; step < 25; ++step) {
        solver.step(solver.stableTimeStep());
    }

    std::vector<double> u, v, p;
    for (std::size_t j = 0; j < 6; ++j) {
        for (std::size_t i = 0; i < 8; ++i) {
            u.push_back(solver.u(i, j));
            v.push_back(solver.v(i, j));
            p.push_back(solver.pressure(i, j));
        }
    }

    FieldArchive archive;
    archive.setMetadata("nx", 8.0);
    archive.setMetadata("ny", 6.0);
    archive.setMetadata("lengthX", 1.0);
    archive.setMetadata("lengthY", 0.75);
    archive.setMetadata("viscosity", 0.01);
    archive.setMetadata("lidVelocity", 1.0);
    archive.setMetadata("time", solver.time());
    archive.setField("u", u);
    archive.setField("v", v);
    archive.setField("p", p);

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_persistence_test_cavity.aecf";
    archive.save(path.string());
    const FieldArchive loaded = FieldArchive::load(path.string());
    std::filesystem::remove(path);

    AETHER_CHECK(loaded.metadata("time") == solver.time());
    AETHER_CHECK(loaded.field("u") == u);
    AETHER_CHECK(loaded.field("v") == v);
    AETHER_CHECK(loaded.field("p") == p);

    // Sanity: after 25 steps of a moving lid, the flow must not still be at
    // rest (guards against this test accidentally validating an all-zeros
    // no-op field, the same discipline used in aether_gpu_tests).
    bool anyNonZero = false;
    for (double value : u) {
        if (value != 0.0) {
            anyNonZero = true;
            break;
        }
    }
    AETHER_CHECK(anyNonZero);
}

// Closes the gap testRoundTripsRealSolverState() deliberately left open:
// not just "the archive round-trips the data faithfully", but "resuming a
// live solver from a saved checkpoint reproduces the exact same trajectory
// as never having stopped it". Both branches run the identical
// deterministic arithmetic (stableTimeStep() derives dt from the current
// field state, so a bit-exact restored state must select the same dt and
// therefore take the same step) -- so equality is expected to be exact,
// not approximate, the same reasoning as testGenericRoundTripIsExact()'s
// exact comparison.
void testResumeMatchesUninterruptedRun2D() {
    const std::size_t nx = 8;
    const std::size_t ny = 6;

    LidDrivenCavitySolver2D uninterrupted(nx, ny, 1.0, 0.75, 0.01, 1.0);
    for (int step = 0; step < 60; ++step) {
        uninterrupted.step(uninterrupted.stableTimeStep());
    }

    LidDrivenCavitySolver2D firstHalf(nx, ny, 1.0, 0.75, 0.01, 1.0);
    for (int step = 0; step < 30; ++step) {
        firstHalf.step(firstHalf.stableTimeStep());
    }

    std::vector<double> u, v, p;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            u.push_back(firstHalf.u(i, j));
            v.push_back(firstHalf.v(i, j));
            p.push_back(firstHalf.pressure(i, j));
        }
    }
    FieldArchive archive;
    archive.setMetadata("time", firstHalf.time());
    archive.setField("u", u);
    archive.setField("v", v);
    archive.setField("p", p);

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_persistence_test_resume2d.aecf";
    archive.save(path.string());
    const FieldArchive loaded = FieldArchive::load(path.string());
    std::filesystem::remove(path);

    LidDrivenCavitySolver2D resumed(nx, ny, 1.0, 0.75, 0.01, 1.0);
    resumed.loadState(loaded.field("u"), loaded.field("v"), loaded.field("p"), loaded.metadata("time"));
    for (int step = 0; step < 30; ++step) {
        resumed.step(resumed.stableTimeStep());
    }

    AETHER_CHECK(resumed.time() == uninterrupted.time());
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            AETHER_CHECK(resumed.u(i, j) == uninterrupted.u(i, j));
            AETHER_CHECK(resumed.v(i, j) == uninterrupted.v(i, j));
            AETHER_CHECK(resumed.pressure(i, j) == uninterrupted.pressure(i, j));
        }
    }
}

// Same proof, for the 3D staggered laminar solver -- exercises
// StaggeredCavityBase3D::loadState() (via StaggeredLidDrivenCavitySolver3D's
// `using` declaration), not just LidDrivenCavitySolver2D's own copy.
void testResumeMatchesUninterruptedRun3D() {
    const std::size_t nx = 4;
    const std::size_t ny = 4;
    const std::size_t nz = 3;

    StaggeredLidDrivenCavitySolver3D uninterrupted(nx, ny, nz, 1.0, 1.0, 0.75, 0.01, 1.0);
    for (int step = 0; step < 20; ++step) {
        uninterrupted.step(uninterrupted.stableTimeStep());
    }

    StaggeredLidDrivenCavitySolver3D firstHalf(nx, ny, nz, 1.0, 1.0, 0.75, 0.01, 1.0);
    for (int step = 0; step < 10; ++step) {
        firstHalf.step(firstHalf.stableTimeStep());
    }

    std::vector<double> u, v, w, p;
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i <= nx; ++i) {
                u.push_back(firstHalf.u(i, j, k));
            }
        }
    }
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j <= ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                v.push_back(firstHalf.v(i, j, k));
            }
        }
    }
    for (std::size_t k = 0; k <= nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                w.push_back(firstHalf.w(i, j, k));
            }
        }
    }
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                p.push_back(firstHalf.pressure(i, j, k));
            }
        }
    }

    FieldArchive archive;
    archive.setMetadata("time", firstHalf.time());
    archive.setField("u", u);
    archive.setField("v", v);
    archive.setField("w", w);
    archive.setField("p", p);

    const std::filesystem::path path =
        std::filesystem::temp_directory_path() / "aether_persistence_test_resume3d.aecf";
    archive.save(path.string());
    const FieldArchive loaded = FieldArchive::load(path.string());
    std::filesystem::remove(path);

    StaggeredLidDrivenCavitySolver3D resumed(nx, ny, nz, 1.0, 1.0, 0.75, 0.01, 1.0);
    resumed.loadState(loaded.field("u"), loaded.field("v"), loaded.field("w"), loaded.field("p"),
                       loaded.metadata("time"));
    for (int step = 0; step < 10; ++step) {
        resumed.step(resumed.stableTimeStep());
    }

    AETHER_CHECK(resumed.time() == uninterrupted.time());
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i < nx; ++i) {
                AETHER_CHECK(resumed.pressure(i, j, k) == uninterrupted.pressure(i, j, k));
            }
        }
    }
    for (std::size_t k = 0; k < nz; ++k) {
        for (std::size_t j = 0; j < ny; ++j) {
            for (std::size_t i = 0; i <= nx; ++i) {
                AETHER_CHECK(resumed.u(i, j, k) == uninterrupted.u(i, j, k));
            }
        }
    }
}

// StructuredGrid3D <-> FieldArchive: min()/max() are stored directly
// (rather than re-deriving max from the loaded spacing, which would not be
// guaranteed to invert the original division bit-for-bit -- see
// GridArchive.hpp), so nx/ny/nz, min, max and the derived spacing must all
// come back exactly.
void testGridArchiveRoundTripIsExact() {
    const StructuredGrid3D original(Vector3(-1.0, 0.5, 2.0), Vector3(3.0, 4.5, 7.0), 5, 6, 7);

    FieldArchive archive;
    aether::persistence::saveGrid(archive, original);
    const StructuredGrid3D loaded = aether::persistence::loadGrid(archive);

    AETHER_CHECK(loaded.nx() == original.nx());
    AETHER_CHECK(loaded.ny() == original.ny());
    AETHER_CHECK(loaded.nz() == original.nz());
    AETHER_CHECK(loaded.min().x == original.min().x && loaded.min().y == original.min().y &&
                 loaded.min().z == original.min().z);
    AETHER_CHECK(loaded.max().x == original.max().x && loaded.max().y == original.max().y &&
                 loaded.max().z == original.max().z);
    AETHER_CHECK(loaded.spacing().x == original.spacing().x && loaded.spacing().y == original.spacing().y &&
                 loaded.spacing().z == original.spacing().z);
}

// ProjectHistory: three checkpoints recorded through one instance must be
// visible, in order, with their saved field data intact, through a
// *second*, independently-constructed instance pointed at the same
// directory -- proving the index is genuinely persisted to disk, not just
// held in memory (a bug here would only be caught by re-opening, which a
// same-instance test can't exercise).
void testProjectHistoryPersistsAcrossInstances() {
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "aether_persistence_test_history";
    std::filesystem::remove_all(dir);

    {
        ProjectHistory history(dir.string());
        for (int i = 0; i < 3; ++i) {
            FieldArchive archive;
            archive.setMetadata("index", static_cast<double>(i));
            archive.setField("u", {static_cast<double>(i), static_cast<double>(i) + 0.5});
            history.record("run_" + std::to_string(i), archive);
        }
    }

    ProjectHistory reopened(dir.string());
    const std::vector<HistoryEntry> entries = reopened.entries();
    AETHER_CHECK(entries.size() == 3);
    for (std::size_t i = 0; i < entries.size(); ++i) {
        AETHER_CHECK(entries[i].label == "run_" + std::to_string(i));
        const FieldArchive archive = reopened.load(entries[i]);
        AETHER_CHECK(archive.metadata("index") == static_cast<double>(i));
        const std::vector<double> expected = {static_cast<double>(i), static_cast<double>(i) + 0.5};
        AETHER_CHECK(archive.field("u") == expected);
    }

    std::filesystem::remove_all(dir);
}

} // namespace

int main() {
    testGenericRoundTripIsExact();
    testLoadRejectsBadMagic();
    testMissingKeysThrow();
    testRoundTripsRealSolverState();
    testResumeMatchesUninterruptedRun2D();
    testResumeMatchesUninterruptedRun3D();
    testGridArchiveRoundTripIsExact();
    testProjectHistoryPersistsAcrossInstances();
    std::printf("aether_persistence_tests: OK\n");
    return 0;
}
