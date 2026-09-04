#pragma once

#include "aether/geometry/TriangleMesh.hpp"

#include <string>
#include <vector>

namespace aether::geometry {

// STEP (ISO 10303-21) import.
//
// **Two build configurations, one contract.** `loadStep()`'s signature and
// `StepLoadResult` never change -- what changes is which implementation
// answers them, selected entirely at compile time (see StepIO.cpp; the
// macro that picks between them, AETHER_HAVE_OPENCASCADE, never appears in
// this header, the same ABI-stability rule engine/solver's own
// AETHER_HAVE_GPU already follows, so no caller of this header can ever
// observe a difference in this class's *shape*, only in its *capability*):
//
//   - **Without OpenCASCADE** (no external dependency): **the
//     planar-faceted solid representation only**. A STEP file's geometry
//     divides cleanly into two kinds, verified against ISO 10303-42's own
//     published EXPRESS schema (steptools.com's hosted SMRL documentation)
//     before writing a line of extraction code here, rather than assumed
//     from memory the way this project avoids everywhere else:
//
//       - A **faceted boundary representation**: every face is planar and
//         every loop is a POLY_LOOP -- literally a flat ordered list of
//         CARTESIAN_POINTs, with no curve or surface entity involved at
//         all. FACETED_BREP is defined in the standard as exactly the
//         restriction of MANIFOLD_SOLID_BREP to this case. This is already
//         triangulable geometry; loading it needs a text-format parser and
//         a planar triangulator, both of which this engine already has.
//
//       - Everything else -- a curved surface (cylindrical, B-spline,
//         toroidal, ...), a trimmed/advanced face, or a tessellated
//         (AP242) shape representation -- is reported as unsupported
//         rather than silently producing an empty or partial mesh: see
//         `StepLoadResult::unsupportedFeatures`.
//
//     Verified during development, not assumed: POLY_LOOP, FACE_BOUND /
//     FACE_OUTER_BOUND, FACE, CLOSED_SHELL and (FACETED_)MANIFOLD_SOLID_BREP's
//     exact attribute lists and order were each checked against the
//     published schema. The one representation that was *not* verified
//     closely enough to implement with the same confidence -- AP242's
//     tessellated TRIANGULATED_SURFACE_SET, which STEP instantiates via
//     Part 21's complex-entity syntax rather than a single simple entity --
//     is deliberately left as a detected-but-unsupported case instead of a
//     guessed implementation. A face bound with more than one loop (an
//     inner loop, i.e. a hole) is also reported as unsupported:
//     PolygonTriangulation2D, reused here for each planar face, does not
//     yet triangulate polygons with holes.
//
//   - **With OpenCASCADE** (ROADMAP Fase 5, `AETHER_ENABLE_OPENCASCADE` in
//     the root CMakeLists.txt, detected via vcpkg the same "opt-in, not
//     required" way `AETHER_ENABLE_CUDA` already is): a real CAD kernel
//     reads the *entire* file -- curved and trimmed faces included -- and
//     tessellates it (`BRepMesh_IncrementalMesh`). There is no longer a
//     category of "valid geometry this loader cannot interpret", so a file
//     the kernel genuinely cannot resolve throws rather than populating
//     `unsupportedFeatures`. `loadIges()` (below) only exists in this
//     configuration in spirit -- it is declared unconditionally, for the
//     ABI reason above, but throws `std::runtime_error` when called in a
//     build without OpenCASCADE, since there is no partial/faceted-only
//     fallback for IGES the way there is for STEP.
struct StepLoadResult {
    TriangleMesh mesh;

    // Non-empty means `mesh` may be incomplete: some entity this file
    // referenced could not be interpreted by the rules above. Each entry
    // names the STEP entity id and the reason, e.g. "#42: B_SPLINE_SURFACE
    // (curved face -- needs a CAD kernel)". An empty mesh with a non-empty
    // list means the file's solid is not faceted at all; a partial mesh
    // with entries means some faces loaded and others did not.
    std::vector<std::string> unsupportedFeatures;
};

// True when this library was built with OpenCASCADE support (Fase 5) --
// the runtime-queryable form of AETHER_HAVE_OPENCASCADE for a caller (a
// Python test, above all, which has no access to a compile-time macro)
// that needs to know which of loadStep()'s two implementations is active
// before choosing what to expect from it. Always declared, same ABI
// reasoning as loadIges() below and as StaggeredCavityBase3D::gpuActive()
// elsewhere in this project.
bool stepIoHasOpenCascade();

// Parses `path` as an ISO 10303-21 (Part 21) file and extracts every
// FACETED_BREP / MANIFOLD_SOLID_BREP it can fully interpret under the
// restriction documented above.
//
// Vertices are shared from the STEP entity graph itself (each
// CARTESIAN_POINT becomes exactly one TriangleMesh vertex, reused by every
// face that references it), so unlike `loadStl` this does not need
// `weldVertices()` afterwards -- the same reasoning `loadObj` already
// documents for its own already-indexed input.
//
// Throws std::runtime_error if the file cannot be opened, is not
// recognisable as an ISO-10303-21 file (missing the `ISO-10303-21;`
// header), or is malformed Part 21 syntax (unbalanced parentheses or
// quotes, a reference to an entity id that was never defined). A
// structurally valid file whose *geometry* this loader cannot interpret is
// not an error -- that is what `unsupportedFeatures` is for.
StepLoadResult loadStep(const std::string& path);

// Parses `path` as an IGES file via OpenCASCADE (ROADMAP Fase 5). Always
// declared regardless of build configuration -- see this header's own
// class comment above for why -- but throws std::runtime_error in a build
// without OpenCASCADE, since IGES has no faceted-only fallback the way
// loadStep() does.
StepLoadResult loadIges(const std::string& path);

} // namespace aether::geometry
