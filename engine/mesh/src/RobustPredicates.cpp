#include "aether/mesh/RobustPredicates.hpp"

#include "aether/core/BigInt.hpp"

#include <array>
#include <cmath>
#include <cstdint>

namespace aether::mesh {

using aether::core::BigInt;
using aether::core::Vector3;

namespace {

// An exact dyadic rational: value == significand * 2^exponent.
struct Dyadic {
    BigInt significand;
    int exponent = 0;
};

// Every double is *exactly* significand * 2^exponent for some integer
// significand with at most 53 bits and some integer exponent -- this
// conversion is exact, not an approximation (std::frexp/std::ldexp never
// round here: m*2^53 is already an integer for any double m in [0.5,1)).
Dyadic exact(double d) {
    if (d == 0.0) {
        return Dyadic{BigInt::fromInt64(0), 0};
    }
    int exp = 0;
    const double m = std::frexp(d, &exp); // d = m * 2^exp, 0.5 <= |m| < 1
    const double scaled = std::ldexp(m, 53); // exact: an integer value
    const int64_t mantissa = static_cast<int64_t>(scaled);
    return Dyadic{BigInt::fromInt64(mantissa), exp - 53};
}

Dyadic add(const Dyadic& a, const Dyadic& b) {
    const int minExp = a.exponent < b.exponent ? a.exponent : b.exponent;
    const BigInt aShifted = a.significand.shiftLeft(a.exponent - minExp);
    const BigInt bShifted = b.significand.shiftLeft(b.exponent - minExp);
    return Dyadic{aShifted + bShifted, minExp};
}

Dyadic sub(const Dyadic& a, const Dyadic& b) { return add(a, Dyadic{-b.significand, b.exponent}); }

Dyadic mul(const Dyadic& a, const Dyadic& b) { return Dyadic{a.significand * b.significand, a.exponent + b.exponent}; }

int signOf(const Dyadic& d) { return d.significand.sign(); }

using Vec3D = std::array<Dyadic, 3>;

Vec3D exactVec(const Vector3& v) { return {exact(v.x), exact(v.y), exact(v.z)}; }

Vec3D subVec(const Vec3D& a, const Vec3D& b) { return {sub(a[0], b[0]), sub(a[1], b[1]), sub(a[2], b[2])}; }

Dyadic dot(const Vec3D& a, const Vec3D& b) {
    return add(add(mul(a[0], b[0]), mul(a[1], b[1])), mul(a[2], b[2]));
}

Vec3D cross(const Vec3D& a, const Vec3D& b) {
    return {sub(mul(a[1], b[2]), mul(a[2], b[1])), sub(mul(a[2], b[0]), mul(a[0], b[2])),
            sub(mul(a[0], b[1]), mul(a[1], b[0]))};
}

Dyadic normSquared(const Vec3D& v) { return dot(v, v); }

// Same cofactor expansion as the plain-double det3x3() this mirrors:
// r0[0]*(r1[1]*r2[2]-r1[2]*r2[1]) - r0[1]*(r1[0]*r2[2]-r1[2]*r2[0]) +
// r0[2]*(r1[0]*r2[1]-r1[1]*r2[0]).
Dyadic det3(const Vec3D& r0, const Vec3D& r1, const Vec3D& r2) {
    const Dyadic t0 = mul(r0[0], sub(mul(r1[1], r2[2]), mul(r1[2], r2[1])));
    const Dyadic t1 = mul(r0[1], sub(mul(r1[0], r2[2]), mul(r1[2], r2[0])));
    const Dyadic t2 = mul(r0[2], sub(mul(r1[0], r2[1]), mul(r1[1], r2[0])));
    return sub(add(t0, t2), t1);
}

} // namespace

int orientation3D(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    const Vec3D ea = exactVec(a);
    const Vec3D u = subVec(exactVec(b), ea);
    const Vec3D v = subVec(exactVec(c), ea);
    const Vec3D w = subVec(exactVec(d), ea);
    return signOf(dot(u, cross(v, w)));
}

int inSphere3D(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d, const Vector3& p) {
    const Vec3D ep = exactVec(p);
    const Vec3D A = subVec(exactVec(a), ep);
    const Vec3D B = subVec(exactVec(b), ep);
    const Vec3D C = subVec(exactVec(c), ep);
    const Vec3D D = subVec(exactVec(d), ep);

    const Dyadic aNorm = normSquared(A);
    const Dyadic bNorm = normSquared(B);
    const Dyadic cNorm = normSquared(C);
    const Dyadic dNorm = normSquared(D);

    const Dyadic detBCD = det3(B, C, D);
    const Dyadic detACD = det3(A, C, D);
    const Dyadic detABD = det3(A, B, D);
    const Dyadic detABC = det3(A, B, C);

    const Dyadic result =
        sub(add(sub(mul(aNorm, detBCD), mul(bNorm, detACD)), mul(cNorm, detABD)), mul(dNorm, detABC));
    return signOf(result);
}

int orientation2D(const Vector3& a, const Vector3& b, const Vector3& c) {
    const Dyadic ax = exact(a.x);
    const Dyadic ay = exact(a.y);
    const Dyadic bax = sub(exact(b.x), ax);
    const Dyadic bay = sub(exact(b.y), ay);
    const Dyadic cax = sub(exact(c.x), ax);
    const Dyadic cay = sub(exact(c.y), ay);
    return signOf(sub(mul(bax, cay), mul(bay, cax)));
}

int inCircle2D(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& p) {
    const Dyadic px = exact(p.x);
    const Dyadic py = exact(p.y);
    const Dyadic ax = sub(exact(a.x), px);
    const Dyadic ay = sub(exact(a.y), py);
    const Dyadic bx = sub(exact(b.x), px);
    const Dyadic by = sub(exact(b.y), py);
    const Dyadic cx = sub(exact(c.x), px);
    const Dyadic cy = sub(exact(c.y), py);

    const Dyadic aNorm = add(mul(ax, ax), mul(ay, ay));
    const Dyadic bNorm = add(mul(bx, bx), mul(by, by));
    const Dyadic cNorm = add(mul(cx, cx), mul(cy, cy));

    // Same cofactor expansion as the plain-double version this mirrors:
    // aNorm*(bx*cy - cx*by) - bNorm*(ax*cy - cx*ay) + cNorm*(ax*by - bx*ay).
    const Dyadic t0 = mul(aNorm, sub(mul(bx, cy), mul(cx, by)));
    const Dyadic t1 = mul(bNorm, sub(mul(ax, cy), mul(cx, ay)));
    const Dyadic t2 = mul(cNorm, sub(mul(ax, by), mul(bx, ay)));
    return signOf(sub(add(t0, t2), t1));
}

} // namespace aether::mesh
