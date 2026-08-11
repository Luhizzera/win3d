#pragma once

#include "aether/geometry/TriangleMesh.hpp"

#include <string>

namespace aether::geometry {

// Loads an STL file (auto-detects ASCII vs binary) and welds coincident
// vertices (see TriangleMesh::weldVertices) so the result has real shared
// connectivity, not the disconnected per-triangle vertices STL stores on
// disk.
TriangleMesh loadStl(const std::string& path, double weldTolerance = 1e-6);

// Writes a binary STL. Face normals are recomputed from vertex winding, not
// read from the mesh, since TriangleMesh does not store them separately.
void saveStlBinary(const TriangleMesh& mesh, const std::string& path);

} // namespace aether::geometry
