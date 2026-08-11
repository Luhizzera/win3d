#pragma once

#include "aether/geometry/TriangleMesh.hpp"

#include <string>

namespace aether::geometry {

// Loads a Wavefront OBJ file: `v x y z` vertex lines and `f ...` face lines
// (each face token may be a bare vertex index, `v/vt`, `v/vt/vn`, or
// `v//vn` -- only the leading vertex index is used, texture/normal indices
// are read past and discarded). Faces with more than 3 vertices are fan-
// triangulated from the first vertex (v0,v1,v2), (v0,v2,v3), ... -- correct
// for the common case of planar convex faces (e.g. quads from most
// exporters), not a general polygon triangulator for arbitrary non-convex
// ngons; PolygonTriangulation2D exists for that but operates on a 2D
// boundary, not an arbitrarily-oriented 3D planar face, so folding it in
// here would need a projection step this first pass doesn't attempt.
// Only positive (absolute) vertex indices are supported; OBJ's relative
// negative-index form is not handled.
//
// Unlike loadStl(), this does *not* call weldVertices(): OBJ's `v`/`f`
// format already stores one shared vertex list referenced by index, not
// STL's independent per-triangle copies, so there is nothing to weld
// (fixing genuinely duplicated positions that a file happens to list
// twice is a data-quality concern outside this loader's scope, same as
// StlIO not fixing malformed STL).
//
// Vertex normals, texture coordinates, groups/objects, and materials
// (vn/vt/o/g/mtllib/usemtl) are read past and ignored, matching the STL
// loader's precedent of not carrying anything beyond raw geometry.
TriangleMesh loadObj(const std::string& path);

// Writes vertex (`v`) and triangle face (`f`) lines only -- no normals or
// texture coordinates, since TriangleMesh does not store them.
void saveObj(const TriangleMesh& mesh, const std::string& path);

} // namespace aether::geometry
