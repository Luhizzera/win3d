#include "aether/geometry/ObjIO.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace aether::geometry {

using aether::core::Vector3;

namespace {

// Parses the leading vertex index out of a face token ("N", "N/N",
// "N/N/N", or "N//N"), returning it 0-based. Throws on a negative
// (relative) index -- not supported, see the header comment.
std::size_t parseFaceVertexIndex(const std::string& token) {
    const auto slashPos = token.find('/');
    const std::string indexStr = slashPos == std::string::npos ? token : token.substr(0, slashPos);
    const long long oneBasedIndex = std::stoll(indexStr);
    if (oneBasedIndex <= 0) {
        throw std::runtime_error("loadObj: relative/negative face indices are not supported: " + token);
    }
    return static_cast<std::size_t>(oneBasedIndex - 1);
}

} // namespace

TriangleMesh loadObj(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("loadObj: could not open file: " + path);
    }

    TriangleMesh mesh;
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::string tag;
        lineStream >> tag;

        if (tag == "v") {
            double x = 0.0, y = 0.0, z = 0.0;
            lineStream >> x >> y >> z;
            mesh.addVertex(Vector3(x, y, z));
        } else if (tag == "f") {
            std::vector<std::size_t> faceVertices;
            std::string token;
            while (lineStream >> token) {
                faceVertices.push_back(parseFaceVertexIndex(token));
            }
            // Fan-triangulate faces with more than 3 vertices (see header
            // comment for the convex-face caveat).
            for (std::size_t i = 1; i + 1 < faceVertices.size(); ++i) {
                mesh.addTriangle(faceVertices[0], faceVertices[i], faceVertices[i + 1]);
            }
        }
        // Everything else (vn, vt, o, g, mtllib, usemtl, s, comments) is
        // intentionally ignored.
    }
    return mesh;
}

void saveObj(const TriangleMesh& mesh, const std::string& path) {
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("saveObj: could not open file for writing: " + path);
    }

    file << "# Aether CFD Engine - OBJ export\n";
    for (std::size_t i = 0; i < mesh.vertexCount(); ++i) {
        const Vector3& v = mesh.vertex(i);
        file << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }
    for (std::size_t i = 0; i < mesh.triangleCount(); ++i) {
        const Triangle& t = mesh.triangle(i);
        // OBJ indices are 1-based.
        file << "f " << (t.vertices[0] + 1) << " " << (t.vertices[1] + 1) << " " << (t.vertices[2] + 1)
             << "\n";
    }
}

} // namespace aether::geometry
