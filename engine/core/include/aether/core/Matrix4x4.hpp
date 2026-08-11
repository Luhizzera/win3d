#pragma once

#include "aether/core/Vector3.hpp"

#include <array>
#include <cmath>

namespace aether::core {

// Column-major 4x4 matrix (OpenGL/GLSL convention: m[col*4+row]), suitable
// for glUniformMatrix4fv(location, 1, GL_FALSE, matrix.data()) with no
// transposition. Added for Module 8's modern (shader-based) rendering
// pipeline, which has no GLU (gluPerspective/gluLookAt do not exist in a
// core-profile OpenGL context) to build projection/view matrices with.
class Matrix4x4 {
public:
    std::array<float, 16> m{};

    Matrix4x4() { m.fill(0.0f); }

    static Matrix4x4 identity() {
        Matrix4x4 r;
        r.m[0] = 1.0f;
        r.m[5] = 1.0f;
        r.m[10] = 1.0f;
        r.m[15] = 1.0f;
        return r;
    }

    // Right-handed perspective projection matching gluPerspective's
    // convention (OpenGL NDC z in [-1, 1]).
    static Matrix4x4 perspective(double fovYRadians, double aspect, double nearZ, double farZ) {
        Matrix4x4 r;
        const double f = 1.0 / std::tan(fovYRadians / 2.0);
        r.m[0] = static_cast<float>(f / aspect);
        r.m[5] = static_cast<float>(f);
        r.m[10] = static_cast<float>((farZ + nearZ) / (nearZ - farZ));
        r.m[11] = -1.0f;
        r.m[14] = static_cast<float>((2.0 * farZ * nearZ) / (nearZ - farZ));
        return r;
    }

    // Orthographic projection matching glOrtho's convention.
    static Matrix4x4 ortho(double left, double right, double bottom, double top, double nearZ,
                            double farZ) {
        Matrix4x4 r = identity();
        r.m[0] = static_cast<float>(2.0 / (right - left));
        r.m[5] = static_cast<float>(2.0 / (top - bottom));
        r.m[10] = static_cast<float>(-2.0 / (farZ - nearZ));
        r.m[12] = static_cast<float>(-(right + left) / (right - left));
        r.m[13] = static_cast<float>(-(top + bottom) / (top - bottom));
        r.m[14] = static_cast<float>(-(farZ + nearZ) / (farZ - nearZ));
        return r;
    }

    // Right-handed view matrix matching gluLookAt's convention.
    static Matrix4x4 lookAt(const Vector3& eye, const Vector3& center, const Vector3& up) {
        const Vector3 f = (center - eye).normalized();
        const Vector3 s = f.cross(up).normalized();
        const Vector3 u = s.cross(f);

        Matrix4x4 r = identity();
        r.m[0] = static_cast<float>(s.x);
        r.m[4] = static_cast<float>(s.y);
        r.m[8] = static_cast<float>(s.z);
        r.m[1] = static_cast<float>(u.x);
        r.m[5] = static_cast<float>(u.y);
        r.m[9] = static_cast<float>(u.z);
        r.m[2] = static_cast<float>(-f.x);
        r.m[6] = static_cast<float>(-f.y);
        r.m[10] = static_cast<float>(-f.z);
        r.m[12] = static_cast<float>(-s.dot(eye));
        r.m[13] = static_cast<float>(-u.dot(eye));
        r.m[14] = static_cast<float>(f.dot(eye));
        return r;
    }

    // Standard column-major matrix product: (A*B) applies B's transform
    // first, then A's -- e.g. mvp = projection * view * model.
    Matrix4x4 operator*(const Matrix4x4& o) const {
        Matrix4x4 r;
        for (int col = 0; col < 4; ++col) {
            for (int row = 0; row < 4; ++row) {
                float sum = 0.0f;
                for (int k = 0; k < 4; ++k) {
                    sum += m[static_cast<std::size_t>(k * 4 + row)] * o.m[static_cast<std::size_t>(col * 4 + k)];
                }
                r.m[static_cast<std::size_t>(col * 4 + row)] = sum;
            }
        }
        return r;
    }

    const float* data() const { return m.data(); }
};

// Transforms a point (implicit w=1) and perspective-divides by the
// resulting w -- mainly useful for testing a projection/view matrix
// without a GPU round-trip.
inline Vector3 transformPoint(const Matrix4x4& mat, const Vector3& p) {
    const auto& m = mat.m;
    const double x = m[0] * p.x + m[4] * p.y + m[8] * p.z + m[12];
    const double y = m[1] * p.x + m[5] * p.y + m[9] * p.z + m[13];
    const double z = m[2] * p.x + m[6] * p.y + m[10] * p.z + m[14];
    const double w = m[3] * p.x + m[7] * p.y + m[11] * p.z + m[15];
    return Vector3(x / w, y / w, z / w);
}

} // namespace aether::core
