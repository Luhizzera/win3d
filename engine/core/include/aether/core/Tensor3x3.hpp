#pragma once

#include "aether/core/Vector3.hpp"

#include <array>

namespace aether::core {

// Dense 3x3 rank-2 tensor (row-major), used for stress/strain/gradient
// quantities in the physics layer.
class Tensor3x3 {
public:
    Tensor3x3() : m_{{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}} {}

    Tensor3x3(double m00, double m01, double m02,
              double m10, double m11, double m12,
              double m20, double m21, double m22)
        : m_{{{m00, m01, m02}, {m10, m11, m12}, {m20, m21, m22}}} {}

    static Tensor3x3 identity() { return Tensor3x3(1, 0, 0, 0, 1, 0, 0, 0, 1); }
    static Tensor3x3 zero() { return Tensor3x3{}; }

    double& operator()(int row, int col) { return m_[row][col]; }
    double operator()(int row, int col) const { return m_[row][col]; }

    Tensor3x3 operator+(const Tensor3x3& o) const;
    Tensor3x3 operator-(const Tensor3x3& o) const;
    Tensor3x3 operator*(double s) const;
    Tensor3x3 operator*(const Tensor3x3& o) const;
    Vector3 operator*(const Vector3& v) const;

    Tensor3x3 transposed() const;
    double trace() const { return m_[0][0] + m_[1][1] + m_[2][2]; }

private:
    std::array<std::array<double, 3>, 3> m_;
};

} // namespace aether::core
