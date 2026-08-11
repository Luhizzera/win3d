#pragma once

#include <array>
#include <cstdint>
#include <cstddef>

namespace aether::core {

// A fixed-width (1280-bit) signed integer, purpose-built for the exact
// geometric predicates in engine/mesh/RobustPredicates -- NOT a general-
// purpose arbitrary-precision type. 1280 bits is a deliberately generous
// fixed budget: the degree-5 polynomials those predicates evaluate (see
// RobustPredicates.hpp) need at most ~270 bits for realistic double
// inputs (5 factors of up to 53 mantissa bits each), so this leaves
// roughly 4x headroom rather than being tuned to the exact minimum.
//
// Limbs are uint32_t (not uint64_t) specifically so multiplication never
// needs a 128-bit intermediate type -- MSVC has no portable __int128, and
// this project avoids compiler-specific intrinsics where a portable
// alternative (every uint32*uint32 product fits exactly in a uint64_t) is
// just as simple. Stored little-endian (limbs[0] is least significant).
//
// Only the operations RobustPredicates actually needs are implemented:
// construction from a small signed integer, add, subtract, negate,
// multiply, shift-left (by a bit count, used to align dyadic-rational
// terms to a common exponent before adding), and sign(). This is
// deliberately not a reusable general bignum library.
class BigInt {
public:
    static constexpr std::size_t kLimbCount = 40; // 1280 bits

    BigInt() = default;
    static BigInt fromInt64(int64_t v);

    BigInt operator+(const BigInt& o) const;
    BigInt operator-(const BigInt& o) const;
    BigInt operator-() const;
    BigInt operator*(const BigInt& o) const;
    BigInt shiftLeft(int bits) const;

    // -1, 0, or +1.
    int sign() const;
    bool isZero() const;

private:
    bool negative_ = false;
    std::array<uint32_t, kLimbCount> limbs_{}; // magnitude, little-endian

    static std::array<uint32_t, kLimbCount> addMagnitudes(const std::array<uint32_t, kLimbCount>& a,
                                                            const std::array<uint32_t, kLimbCount>& b);
    // Requires magnitude(a) >= magnitude(b).
    static std::array<uint32_t, kLimbCount> subMagnitudes(const std::array<uint32_t, kLimbCount>& a,
                                                            const std::array<uint32_t, kLimbCount>& b);
    static std::array<uint32_t, kLimbCount> mulMagnitudes(const std::array<uint32_t, kLimbCount>& a,
                                                            const std::array<uint32_t, kLimbCount>& b);
    static int compareMagnitudes(const std::array<uint32_t, kLimbCount>& a,
                                  const std::array<uint32_t, kLimbCount>& b);
    static bool isZeroMagnitude(const std::array<uint32_t, kLimbCount>& a);
};

} // namespace aether::core
