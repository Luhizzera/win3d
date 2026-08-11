#include "aether/geometry/StlIO.hpp"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <vector>

namespace aether::geometry {

using aether::core::Vector3;

namespace {

std::uint32_t readBinaryTriangleCount(std::ifstream& file) {
    file.seekg(80, std::ios::beg);
    std::uint32_t count = 0;
    file.read(reinterpret_cast<char*>(&count), sizeof(count));
    return count;
}

TriangleMesh loadBinaryStl(std::ifstream& file, std::uint32_t triangleCount) {
    TriangleMesh mesh;
    file.seekg(84, std::ios::beg);
    for (std::uint32_t i = 0; i < triangleCount; ++i) {
        float normal[3];
        file.read(reinterpret_cast<char*>(normal), sizeof(normal));

        std::size_t indices[3];
        for (int v = 0; v < 3; ++v) {
            float p[3];
            file.read(reinterpret_cast<char*>(p), sizeof(p));
            indices[v] = mesh.addVertex(Vector3(p[0], p[1], p[2]));
        }

        std::uint16_t attributeByteCount = 0;
        file.read(reinterpret_cast<char*>(&attributeByteCount), sizeof(attributeByteCount));

        mesh.addTriangle(indices[0], indices[1], indices[2]);
    }
    return mesh;
}

TriangleMesh loadAsciiStl(std::ifstream& file) {
    file.clear();
    file.seekg(0, std::ios::beg);
    TriangleMesh mesh;

    std::string token;
    std::vector<Vector3> pending;
    pending.reserve(3);

    while (file >> token) {
        if (token == "vertex") {
            double x = 0.0, y = 0.0, z = 0.0;
            file >> x >> y >> z;
            pending.emplace_back(x, y, z);
            if (pending.size() == 3) {
                const std::size_t a = mesh.addVertex(pending[0]);
                const std::size_t b = mesh.addVertex(pending[1]);
                const std::size_t c = mesh.addVertex(pending[2]);
                mesh.addTriangle(a, b, c);
                pending.clear();
            }
        }
    }
    return mesh;
}

} // namespace

TriangleMesh loadStl(const std::string& path, double weldTolerance) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("loadStl: could not open file: " + path);
    }

    file.seekg(0, std::ios::end);
    const std::streamoff fileSize = file.tellg();

    TriangleMesh mesh;
    if (fileSize >= 84) {
        const std::uint32_t triangleCount = readBinaryTriangleCount(file);
        const std::streamoff expectedBinarySize =
            84 + static_cast<std::streamoff>(triangleCount) * 50;
        if (expectedBinarySize == fileSize) {
            mesh = loadBinaryStl(file, triangleCount);
        } else {
            mesh = loadAsciiStl(file);
        }
    } else {
        mesh = loadAsciiStl(file);
    }

    mesh.weldVertices(weldTolerance);
    return mesh;
}

void saveStlBinary(const TriangleMesh& mesh, const std::string& path) {
    std::ofstream file(path, std::ios::binary);
    if (!file) {
        throw std::runtime_error("saveStlBinary: could not open file for writing: " + path);
    }

    char header[80] = {};
    std::strncpy(header, "Aether CFD Engine - binary STL", sizeof(header) - 1);
    file.write(header, sizeof(header));

    const std::uint32_t triangleCount = static_cast<std::uint32_t>(mesh.triangleCount());
    file.write(reinterpret_cast<const char*>(&triangleCount), sizeof(triangleCount));

    for (std::size_t i = 0; i < mesh.triangleCount(); ++i) {
        const Vector3 normal = mesh.faceNormal(i);
        const float n[3] = {
            static_cast<float>(normal.x), static_cast<float>(normal.y), static_cast<float>(normal.z)};
        file.write(reinterpret_cast<const char*>(n), sizeof(n));

        const Triangle& t = mesh.triangle(i);
        for (std::size_t v = 0; v < 3; ++v) {
            const Vector3& p = mesh.vertex(t.vertices[v]);
            const float coords[3] = {
                static_cast<float>(p.x), static_cast<float>(p.y), static_cast<float>(p.z)};
            file.write(reinterpret_cast<const char*>(coords), sizeof(coords));
        }

        const std::uint16_t attributeByteCount = 0;
        file.write(reinterpret_cast<const char*>(&attributeByteCount), sizeof(attributeByteCount));
    }
}

} // namespace aether::geometry
