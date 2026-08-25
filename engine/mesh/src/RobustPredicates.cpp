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

// -- Floating-point filters ----------------------------------------------
//
// Each filter evaluates the very same determinant in plain doubles and,
// alongside it, the expression's **permanent**: the identical formula with
// every operand replaced by its magnitude and every subtraction by an
// addition. The permanent is an upper bound on how large the intermediate
// quantities can get, so C * epsilon * permanent bounds the total roundoff
// the double evaluation can have accumulated, for a constant C covering
// the number of rounding operations along the longest path. If the
// computed determinant exceeds that bound in magnitude, no possible
// roundoff could have flipped its sign, and the double answer is provably
// the exact answer's sign.
//
// The constants below are chosen several times larger than the tightest
// error analysis requires (a degree-3 determinant reached through ~9
// roundings, a degree-5 one through ~20), because the failure modes are
// wildly asymmetric: too large merely routes extra cases through exact
// arithmetic, too small returns a wrong sign. Being generous is close to
// free in practice -- on non-degenerate input the determinant is orders of
// magnitude away from the bound either way, so the constant barely moves
// the fallback rate.
//
// kFilterInconclusive is deliberately not +1/0/-1: it means "the double
// path declined to answer", not a geometric result.
constexpr double kEpsilon = 1.1102230246251565e-16; // 2^-53, double unit roundoff
constexpr double kOrientationFilterConstant = 32.0;
constexpr double kInSphereFilterConstant = 128.0;
constexpr int kFilterInconclusive = 2;

int orientation3DFilter(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    const double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
    const double vx = c.x - a.x, vy = c.y - a.y, vz = c.z - a.z;
    const double wx = d.x - a.x, wy = d.y - a.y, wz = d.z - a.z;

    const double m0 = vy * wz - vz * wy;
    const double m1 = vx * wz - vz * wx;
    const double m2 = vx * wy - vy * wx;
    const double det = ux * m0 - uy * m1 + uz * m2;

    const double permanent = std::fabs(ux) * (std::fabs(vy * wz) + std::fabs(vz * wy)) +
                              std::fabs(uy) * (std::fabs(vx * wz) + std::fabs(vz * wx)) +
                              std::fabs(uz) * (std::fabs(vx * wy) + std::fabs(vy * wx));
    const double bound = kOrientationFilterConstant * kEpsilon * permanent;

    // A zero permanent forces the inconclusive path rather than reporting
    // sign 0: it means every term vanished in double arithmetic, which is
    // exactly the degenerate case worth confirming exactly.
    if (det > bound && det > 0.0) {
        return 1;
    }
    if (det < -bound && det < 0.0) {
        return -1;
    }
    return kFilterInconclusive;
}

int inSphere3DFilter(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d,
                      const Vector3& p) {
    const double ax = a.x - p.x, ay = a.y - p.y, az = a.z - p.z;
    const double bx = b.x - p.x, by = b.y - p.y, bz = b.z - p.z;
    const double cx = c.x - p.x, cy = c.y - p.y, cz = c.z - p.z;
    const double dx = d.x - p.x, dy = d.y - p.y, dz = d.z - p.z;

    const double aNorm = ax * ax + ay * ay + az * az;
    const double bNorm = bx * bx + by * by + bz * bz;
    const double cNorm = cx * cx + cy * cy + cz * cz;
    const double dNorm = dx * dx + dy * dy + dz * dz;

    const double bcd = bx * (cy * dz - cz * dy) - by * (cx * dz - cz * dx) + bz * (cx * dy - cy * dx);
    const double acd = ax * (cy * dz - cz * dy) - ay * (cx * dz - cz * dx) + az * (cx * dy - cy * dx);
    const double abd = ax * (by * dz - bz * dy) - ay * (bx * dz - bz * dx) + az * (bx * dy - by * dx);
    const double abc = ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);

    const double det = aNorm * bcd - bNorm * acd + cNorm * abd - dNorm * abc;

    const auto det3Permanent = [](double r0x, double r0y, double r0z, double r1x, double r1y,
                                   double r1z, double r2x, double r2y, double r2z) {
        return std::fabs(r0x) * (std::fabs(r1y * r2z) + std::fabs(r1z * r2y)) +
               std::fabs(r0y) * (std::fabs(r1x * r2z) + std::fabs(r1z * r2x)) +
               std::fabs(r0z) * (std::fabs(r1x * r2y) + std::fabs(r1y * r2x));
    };
    // The four norms are already non-negative, so they need no fabs.
    const double permanent = aNorm * det3Permanent(bx, by, bz, cx, cy, cz, dx, dy, dz) +
                              bNorm * det3Permanent(ax, ay, az, cx, cy, cz, dx, dy, dz) +
                              cNorm * det3Permanent(ax, ay, az, bx, by, bz, dx, dy, dz) +
                              dNorm * det3Permanent(ax, ay, az, bx, by, bz, cx, cy, cz);
    const double bound = kInSphereFilterConstant * kEpsilon * permanent;

    if (det > bound && det > 0.0) {
        return 1;
    }
    if (det < -bound && det < 0.0) {
        return -1;
    }
    return kFilterInconclusive;
}

} // namespace

int orientation3DExact(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    const Vec3D ea = exactVec(a);
    const Vec3D u = subVec(exactVec(b), ea);
    const Vec3D v = subVec(exactVec(c), ea);
    const Vec3D w = subVec(exactVec(d), ea);
    return signOf(dot(u, cross(v, w)));
}

int orientation3D(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d) {
    const int filtered = orientation3DFilter(a, b, c, d);
    if (filtered != kFilterInconclusive) {
        return filtered;
    }
    return orientation3DExact(a, b, c, d);
}

int inSphere3D(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d, const Vector3& p) {
    const int filtered = inSphere3DFilter(a, b, c, d, p);
    if (filtered != kFilterInconclusive) {
        return filtered;
    }
    return inSphere3DExact(a, b, c, d, p);
}

int inSphere3DExact(const Vector3& a, const Vector3& b, const Vector3& c, const Vector3& d, const Vector3& p) {
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
