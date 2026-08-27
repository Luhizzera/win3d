"""Measures whether the limited convection scheme (DIVIDA_TECNICA.md 4.4)
helps the 3D staggered cavity the way it already measurably helps the 2D
one, before changing StaggeredCavityBase3D's default away from Central.

**Why this script exists instead of a straight decision.** The predictor
that all six 3D cavity closures share (StaggeredCavityBase3D) gained
FirstOrderUpwind and LimitedLinearUpwind alongside its original Central
scheme, exposed for measurement on StaggeredLidDrivenCavitySolver3D. Central
stayed the default because every number this base has ever produced must
stay reproducible bit-for-bit until a measurement -- not an assumption --
justifies changing it, the same discipline the 2D port
(LidDrivenCavitySolver2D) already applied.

**Method, mirrored from the 2D measurement.** Cavity at Re=400 (viscosity =
lid_velocity * L / Re), marched to t=12. RMS deviation of u along the
vertical centerline (x = y = L/2, z varying) at two coarse resolutions,
against a shared finer reference -- linearly interpolated onto the
reference's own z grid by hand (no numpy: it is not a dependency this
project otherwise has, and linear interpolation over ~30 points does not
earn adding one), since the resolutions don't share grid points.
2D used n=128 as its reference (cell Reynolds 3.1); an equivalent 3D
resolution is not affordable in reasonable time (cost scales with n^4, not
n^2), so n=32 (cell Reynolds 12.5) is used instead -- coarser than 2D's
reference in absolute terms, but the same relative role: fine enough that
refining further should not flip which scheme is closest to it, which is
the property this comparison actually needs.

Run: python python/research/convection_scheme_3d_cavity.py
"""

import bisect
import sys
import time

sys.path.insert(0, "python")

import aether

Scheme = aether.StaggeredLidDrivenCavitySolver3D.ConvectionScheme
LENGTH = 1.0
LID_VELOCITY = 1.0
REYNOLDS = 400.0
VISCOSITY = LID_VELOCITY * LENGTH / REYNOLDS
MARCH_TO_TIME = 12.0


def centerline_profile(n, scheme):
    solver = aether.StaggeredLidDrivenCavitySolver3D(
        n, n, n, LENGTH, LENGTH, LENGTH, VISCOSITY, LID_VELOCITY, scheme)
    dt = solver.stable_time_step() * 0.5
    steps = int(MARCH_TO_TIME / dt) + 1
    t0 = time.time()
    for _ in range(steps):
        solver.step(dt)
    elapsed = time.time() - t0

    i = n // 2  # u's own x-face grid: exactly x = L/2 for even n
    j = n // 2  # nearest cell-center row to y = L/2
    dz = LENGTH / n
    z = [(k + 0.5) * dz for k in range(n)]
    u = [solver.u(i, j, k) for k in range(n)]
    return z, u, elapsed, solver.max_divergence()


def interp(z, u, query_z):
    """Piecewise-linear interpolation, clamped at the ends -- a stand-in
    for numpy.interp over a handful of points, not a general tool."""
    idx = bisect.bisect_left(z, query_z)
    if idx <= 0:
        return u[0]
    if idx >= len(z):
        return u[-1]
    z0, z1 = z[idx - 1], z[idx]
    u0, u1 = u[idx - 1], u[idx]
    t = (query_z - z0) / (z1 - z0)
    return u0 + t * (u1 - u0)


def rms_against_reference(z, u, ref_z, ref_u):
    squared_errors = [(interp(z, u, q) - r) ** 2 for q, r in zip(ref_z, ref_u)]
    return (sum(squared_errors) / len(squared_errors)) ** 0.5


def main():
    print(f"Re={REYNOLDS:.0f}  viscosity={VISCOSITY:.6f}  marchando a t={MARCH_TO_TIME}")
    print()

    reference_n = 32
    print(f"referencia: n={reference_n} (Re de celula {REYNOLDS / reference_n:.1f}), esquema Central")
    ref_z, ref_u, ref_elapsed, ref_div = centerline_profile(reference_n, Scheme.CENTRAL)
    print(f"  {ref_elapsed:.1f}s, divergencia maxima={ref_div:.3e}")
    print()

    schemes = [
        ("central", Scheme.CENTRAL),
        ("upwind 1a ordem", Scheme.FIRST_ORDER_UPWIND),
        ("limitado (van Leer)", Scheme.LIMITED_LINEAR_UPWIND),
    ]
    test_resolutions = [8, 16]

    header = f"{'esquema':<22}" + "".join(f"n={n:<3d} (ReCel {REYNOLDS / n:5.1f})   " for n in test_resolutions)
    print(header)
    for name, scheme in schemes:
        row = f"{name:<22}"
        for n in test_resolutions:
            z, u, elapsed, div = centerline_profile(n, scheme)
            rms = rms_against_reference(z, u, ref_z, ref_u)
            row += f"rms={rms:.6f} ({elapsed:4.1f}s)  "
            if div > 1e-2:
                row += f"[divergencia alta: {div:.2e}] "
        print(row)


if __name__ == "__main__":
    main()
