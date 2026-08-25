#pragma once

#include "aether/geometry/TriangleMesh.hpp"

#include <string>
#include <vector>

namespace aether::geometry {

// STEP (ISO 10303-21) import -- **the planar-faceted solid representation
// only**, not general CAD.
//
// **Why this scope, and why it is a real boundary rather than a
// placeholder.** A STEP file's geometry divides cleanly into two kinds,
// verified against ISO 10303-42's own published EXPRESS schema (steptools.com's
// hosted SMRL documentation) before writing a line of extraction code here,
// rather than assumed from memory the way this project avoids everywhere
// else:
//
//   - A **faceted boundary representation**: every face is planar and every
//     loop is a POLY_LOOP -- literally a flat ordered list of
//     CARTESIAN_POINTs, with no curve or surface entity involved at all.
//     FACETED_BREP is defined in the standard as exactly the restriction of
//     MANIFOLD_SOLID_BREP to this case. This is already triangulable
//     geometry; loading it needs a text-format parser and a planar
//     triangulator, both of which this engine already has.
//
//   - Everything else -- a curved surface (cylindrical, B-spline, toroidal,
//     ...), a trimmed/advanced face, or a tessellated (AP242) shape
//     representation -- needs an actual CAD kernel to evaluate. This loader
//     does not attempt any of it, and **says so explicitly** rather than
//     silently returning an empty or partial mesh: see
//     `StepLoadResult::unsupportedFeatures`. Building that path is a real,
//     separate decision (OpenCASCADE is the realistic option, and it is a
//     second heavy external dependency after CUDA), left open in
//     ROADMAP.md's own Fase 5 rather than guessed at here.
//
// Verified during development, not assumed: POLY_LOOP, FACE_BOUND /
// FACE_OUTER_BOUND, FACE, CLOSED_SHELL and (FACETED_)MANIFOLD_SOLID_BREP's
// exact attribute lists and order were each checked against the published
// schema. The one representation that was *not* verified closely enough to
// implement with the same confidence -- AP242's tessellated
// TRIANGULATED_SURFACE_SET, which STEP instantiates via Part 21's complex-
// entity syntax rather than a single simple entity -- is deliberately left
// as a detected-but-unsupported case instead of a guessed implementation.
//
// A face bound with more than one loop (an inner loop, i.e. a hole) is
// also reported as unsupported: PolygonTriangulation2D, reused here for
// each planar face, does not yet triangulate polygons with holes.
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

} // namespace aether::geometry
