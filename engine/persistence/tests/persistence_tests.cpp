#include "aether/persistence/FieldArchive.hpp"
#include "aether/persistence/GridArchive.hpp"
#include "aether/persistence/ProjectHistory.hpp"
#include "aether/persistence/TetrahedralMeshArchive.hpp"
#include "aether/solver/LidDrivenCavitySolver2D.hpp"
#include "aether/solver/StaggeredLidDrivenCavitySolver3D.hpp"
#include "aether/testing/Check.hpp"

#include <filesystem>
#include <array>
#include <algorithm>
#include <cmath>
#include <string>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <vector>

using aether::core::Vector3;
using aether::mesh::StructuredGrid3D;
using aether::persistence::FieldArchive;
using aether::mesh::DelaunayTetrahedralization3D;
using aether::mesh::TetrahedralMesh;
using aether::persistence::saveTetrahedralMesh;
using aether::persistence::loadTetrahedralMesh;
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

// **A mesh checkpoint is only worth anything if the reloaded mesh is the
// same mesh**, so this checks identity rather than plausibility: every
// vertex position bit for bit, every cell's connectivity, and -- because
// those two alone would not catch a face-connectivity bug -- the derived
// quantities the finite-volume layer actually consumes: per-cell volume,
// per-cell centroid, face count, and the discrete divergence theorem on
// every cell.
//
// The last of those is the sharpest: a closed polyhedron's outward area
// vectors sum to identically zero, so a reload that lost or misoriented a
// face shows up there even if the vertex and cell arrays round-tripped
// perfectly.
void testTetrahedralMeshArchiveRoundTripsExactly() {
    DelaunayTetrahedralization3D tets;
    // A jittered lattice rather than a regular one, for the reason the mesh
    // suite already documents: a regular lattice tetrahedralizes into
    // co-spherical ties, the degenerate case, which would make this a
    // weaker test than it looks.
    std::uint64_t state = 0x243F6A8885A308D3ull;
    const auto nextJitter = [&state]() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        return (static_cast<double>(state >> 11) / 9007199254740992.0 - 0.5) * 0.2;
    };
    for (int i = 0; i <= 2; ++i) {
        for (int j = 0; j <= 2; ++j) {
            for (int k = 0; k <= 2; ++k) {
                const bool interior = i > 0 && i < 2 && j > 0 && j < 2 && k > 0 && k < 2;
                tets.addPoint(i * 0.5 + (interior ? nextJitter() : 0.0),
                              j * 0.5 + (interior ? nextJitter() : 0.0),
                              k * 0.5 + (interior ? nextJitter() : 0.0));
            }
        }
    }
    tets.tetrahedralize();
    const TetrahedralMesh original = TetrahedralMesh::fromTetrahedralization(tets);
    AETHER_CHECK(original.cellCount() > 0);

    const std::string path = "aether_tet_mesh_roundtrip.aecf";
    {
        FieldArchive archive;
        saveTetrahedralMesh(archive, original);
        // A field defined over the mesh travels in the same archive, which
        // is the whole point of reusing FieldArchive rather than inventing
        // a second container: a checkpoint whose mesh and fields could be
        // separated is a checkpoint that can be silently mismatched.
        std::vector<double> pressure(original.cellCount());
        for (std::size_t c = 0; c < original.cellCount(); ++c) {
            pressure[c] = static_cast<double>(c) * 0.25;
        }
        archive.setField("pressure", pressure);
        archive.save(path);
    }

    const FieldArchive reloaded = FieldArchive::load(path);
    const TetrahedralMesh restored = loadTetrahedralMesh(reloaded);
    std::remove(path.c_str());

    AETHER_CHECK(restored.vertexCount() == original.vertexCount());
    AETHER_CHECK(restored.cellCount() == original.cellCount());
    AETHER_CHECK(restored.faceCount() == original.faceCount());

    for (std::size_t v = 0; v < original.vertexCount(); ++v) {
        // Exact, not near: FieldArchive stores raw doubles, so anything
        // other than bit equality is a bug rather than accumulated error.
        AETHER_CHECK(restored.vertex(v).x == original.vertex(v).x);
        AETHER_CHECK(restored.vertex(v).y == original.vertex(v).y);
        AETHER_CHECK(restored.vertex(v).z == original.vertex(v).z);
    }

    double worstAreaSum = 0.0;
    for (std::size_t c = 0; c < original.cellCount(); ++c) {
        AETHER_CHECK(restored.cellVertices(c) == original.cellVertices(c));
        AETHER_CHECK(restored.cellVolume(c) == original.cellVolume(c));
        AETHER_CHECK(restored.cellCentroid(c).x == original.cellCentroid(c).x);
        AETHER_CHECK(restored.cellCentroid(c).y == original.cellCentroid(c).y);
        AETHER_CHECK(restored.cellCentroid(c).z == original.cellCentroid(c).z);
        worstAreaSum = std::max(worstAreaSum, restored.cellAreaVectorSum(c).norm());
    }
    AETHER_CHECK(worstAreaSum < 1e-12);
    AETHER_CHECK(restored.totalVolume() == original.totalVolume());

    // The co-travelling field survives untouched.
    AETHER_CHECK(reloaded.hasField("pressure"));
    AETHER_CHECK(reloaded.field("pressure").size() == original.cellCount());
    AETHER_CHECK(reloaded.field("pressure")[1] == 0.25);

    std::printf("  [persistence_tests] malha tetraedrica: %zu vertices, %zu celulas, %zu faces "
                "identicas apos ida e volta (pior soma de areas %.2e)\n",
                restored.vertexCount(), restored.cellCount(), restored.faceCount(), worstAreaSum);
    std::fflush(stdout);
}

// A truncated or foreign archive has to fail where it is read. Checked
// because a validation nobody exercises is a validation that quietly stops
// working -- and because the failure it prevents (a mesh silently missing
// its last cell) is exactly the kind that surfaces much later as a wrong
// answer rather than as an error.
void testTetrahedralMeshArchiveRefusesMalformedInput() {
    bool refusedEmpty = false;
    try {
        FieldArchive empty;
        loadTetrahedralMesh(empty);
    } catch (const std::runtime_error&) {
        refusedEmpty = true;
    }
    AETHER_CHECK(refusedEmpty);

    bool refusedTruncated = false;
    try {
        FieldArchive truncated;
        truncated.setField("tetrahedral_mesh_vertices", {0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                                          0.0, 1.0, 0.0, 0.0, 0.0, 1.0});
        truncated.setField("tetrahedral_mesh_cells", {0.0, 1.0, 2.0}); // three, not four
        loadTetrahedralMesh(truncated);
    } catch (const std::runtime_error&) {
        refusedTruncated = true;
    }
    AETHER_CHECK(refusedTruncated);

    bool refusedOutOfRange = false;
    try {
        FieldArchive bad;
        bad.setField("tetrahedral_mesh_vertices", {0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                                                    0.0, 1.0, 0.0, 0.0, 0.0, 1.0});
        bad.setField("tetrahedral_mesh_cells", {0.0, 1.0, 2.0, 99.0}); // no vertex 99
        loadTetrahedralMesh(bad);
    } catch (const std::exception&) {
        refusedOutOfRange = true;
    }
    AETHER_CHECK(refusedOutOfRange);

    std::printf("  [persistence_tests] malha tetraedrica: arquivo vazio, truncado e com indice "
                "fora de faixa recusados\n");
    std::fflush(stdout);
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
    testTetrahedralMeshArchiveRoundTripsExactly();
    testTetrahedralMeshArchiveRefusesMalformedInput();
    std::printf("aether_persistence_tests: OK\n");
    return 0;
}
