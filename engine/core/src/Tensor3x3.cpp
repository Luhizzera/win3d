#include "aether/core/Tensor3x3.hpp"

namespace aether::core {

Tensor3x3 Tensor3x3::operator+(const Tensor3x3& o) const {
    Tensor3x3 result;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            result(r, c) = (*this)(r, c) + o(r, c);
        }
    }
    return result;
}

Tensor3x3 Tensor3x3::operator-(const Tensor3x3& o) const {
    Tensor3x3 result;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            result(r, c) = (*this)(r, c) - o(r, c);
        }
    }
    return result;
}

Tensor3x3 Tensor3x3::operator*(double s) const {
    Tensor3x3 result;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            result(r, c) = (*this)(r, c) * s;
        }
    }
    return result;
}

Tensor3x3 Tensor3x3::operator*(const Tensor3x3& o) const {
    Tensor3x3 result;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k) {
                sum += (*this)(r, k) * o(k, c);
            }
            result(r, c) = sum;
        }
    }
    return result;
}

Vector3 Tensor3x3::operator*(const Vector3& v) const {
    return {
        (*this)(0, 0) * v.x + (*this)(0, 1) * v.y + (*this)(0, 2) * v.z,
        (*this)(1, 0) * v.x + (*this)(1, 1) * v.y + (*this)(1, 2) * v.z,
        (*this)(2, 0) * v.x + (*this)(2, 1) * v.y + (*this)(2, 2) * v.z,
    };
}

Tensor3x3 Tensor3x3::transposed() const {
    Tensor3x3 result;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            result(c, r) = (*this)(r, c);
        }
    }
    return result;
}

} // namespace aether::core
