#include "aether/core/BigInt.hpp"

namespace aether::core {

namespace {
using Limbs = std::array<uint32_t, BigInt::kLimbCount>;
}

BigInt BigInt::fromInt64(int64_t v) {
    BigInt result;
    if (v == 0) {
        return result;
    }
    result.negative_ = v < 0;
    // Convert via uint64_t to sidestep the INT64_MIN negation-overflow
    // corner case (not something this project's mantissa values -- always
    // well within +-2^53 -- ever produce, but correct regardless).
    const uint64_t magnitude = result.negative_ ? (~static_cast<uint64_t>(v) + 1) : static_cast<uint64_t>(v);
    result.limbs_[0] = static_cast<uint32_t>(magnitude & 0xFFFFFFFFu);
    result.limbs_[1] = static_cast<uint32_t>(magnitude >> 32);
    return result;
}

bool BigInt::isZeroMagnitude(const Limbs& a) {
    for (uint32_t limb : a) {
        if (limb != 0) {
            return false;
        }
    }
    return true;
}

int BigInt::compareMagnitudes(const Limbs& a, const Limbs& b) {
    for (std::size_t i = kLimbCount; i-- > 0;) {
        if (a[i] != b[i]) {
            return a[i] < b[i] ? -1 : 1;
        }
    }
    return 0;
}

Limbs BigInt::addMagnitudes(const Limbs& a, const Limbs& b) {
    Limbs result{};
    uint64_t carry = 0;
    for (std::size_t i = 0; i < kLimbCount; ++i) {
        const uint64_t sum = static_cast<uint64_t>(a[i]) + b[i] + carry;
        result[i] = static_cast<uint32_t>(sum & 0xFFFFFFFFu);
        carry = sum >> 32;
    }
    // Any final carry beyond kLimbCount limbs is dropped -- see the header
    // comment: this type's fixed width is sized with generous headroom
    // for its actual use, not a general arbitrary-precision guarantee.
    return result;
}

Limbs BigInt::subMagnitudes(const Limbs& a, const Limbs& b) {
    // Requires magnitude(a) >= magnitude(b) (enforced by the caller).
    Limbs result{};
    int64_t borrow = 0;
    for (std::size_t i = 0; i < kLimbCount; ++i) {
        int64_t diff = static_cast<int64_t>(a[i]) - b[i] - borrow;
        if (diff < 0) {
            diff += (int64_t{1} << 32);
            borrow = 1;
        } else {
            borrow = 0;
        }
        result[i] = static_cast<uint32_t>(diff);
    }
    return result;
}

Limbs BigInt::mulMagnitudes(const Limbs& a, const Limbs& b) {
    // Schoolbook multiplication, truncated to kLimbCount limbs of result.
    // Every partial product a[i]*b[j] fits exactly in a uint64_t (each
    // factor is a uint32_t), so no wider intermediate type is ever needed.
    Limbs result{};
    for (std::size_t i = 0; i < kLimbCount; ++i) {
        if (a[i] == 0) {
            continue;
        }
        uint64_t carry = 0;
        for (std::size_t j = 0; i + j < kLimbCount; ++j) {
            const uint64_t product =
                static_cast<uint64_t>(a[i]) * b[j] + result[i + j] + carry;
            result[i + j] = static_cast<uint32_t>(product & 0xFFFFFFFFu);
            carry = product >> 32;
        }
        // Any carry propagating beyond kLimbCount is dropped (see the
        // header comment on the fixed-width tradeoff).
    }
    return result;
}

BigInt BigInt::operator+(const BigInt& o) const {
    BigInt result;
    if (negative_ == o.negative_) {
        result.negative_ = negative_;
        result.limbs_ = addMagnitudes(limbs_, o.limbs_);
    } else {
        const int cmp = compareMagnitudes(limbs_, o.limbs_);
        if (cmp == 0) {
            return BigInt(); // zero
        }
        if (cmp > 0) {
            result.negative_ = negative_;
            result.limbs_ = subMagnitudes(limbs_, o.limbs_);
        } else {
            result.negative_ = o.negative_;
            result.limbs_ = subMagnitudes(o.limbs_, limbs_);
        }
    }
    if (isZeroMagnitude(result.limbs_)) {
        result.negative_ = false;
    }
    return result;
}

BigInt BigInt::operator-(const BigInt& o) const { return *this + (-o); }

BigInt BigInt::operator-() const {
    BigInt result = *this;
    if (!isZeroMagnitude(result.limbs_)) {
        result.negative_ = !result.negative_;
    }
    return result;
}

BigInt BigInt::operator*(const BigInt& o) const {
    BigInt result;
    result.limbs_ = mulMagnitudes(limbs_, o.limbs_);
    result.negative_ = (negative_ != o.negative_) && !isZeroMagnitude(result.limbs_);
    return result;
}

BigInt BigInt::shiftLeft(int bits) const {
    if (bits <= 0 || isZeroMagnitude(limbs_)) {
        return *this;
    }
    const int limbShift = bits / 32;
    const int bitShift = bits % 32;

    BigInt result;
    result.negative_ = negative_;
    if (static_cast<std::size_t>(limbShift) >= kLimbCount) {
        return result; // shifted entirely out of range -- underflows to 0
    }
    for (std::size_t i = kLimbCount; i-- > static_cast<std::size_t>(limbShift);) {
        const std::size_t src = i - static_cast<std::size_t>(limbShift);
        uint64_t value = static_cast<uint64_t>(limbs_[src]) << bitShift;
        if (bitShift > 0 && src > 0) {
            value |= static_cast<uint64_t>(limbs_[src - 1]) >> (32 - bitShift);
        }
        result.limbs_[i] = static_cast<uint32_t>(value & 0xFFFFFFFFu);
    }
    return result;
}

int BigInt::sign() const {
    if (isZeroMagnitude(limbs_)) {
        return 0;
    }
    return negative_ ? -1 : 1;
}

bool BigInt::isZero() const { return isZeroMagnitude(limbs_); }

} // namespace aether::core
