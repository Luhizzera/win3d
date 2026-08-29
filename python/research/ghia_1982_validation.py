"""Validates LidDrivenCavitySolver2D against Ghia, Ghia & Shin (1982), the
canonical published benchmark for the lid-driven cavity -- external
credibility that no amount of this project's own internal testing can buy
by itself (see the "beta" line of work's item 3).

**Source discipline, the same this project already applies to physics.**
The reference data below was NOT recalled from memory -- transcribing a
34-point table from memory is exactly the risk this project's own
DIVIDA_TECNICA.md names as the reason it avoided literature benchmarks
until now. It was fetched from two independently-maintained transcriptions
of the paper's own Table I (u-velocity) and Table II (v-velocity):

    https://gist.github.com/ivan-pi/3e9326d18a366ffe6a8e5bfda6353219  (u)
    https://gist.github.com/ivan-pi/caa6c6737d36a9140fbcf2ea59c78b3c  (v)

and the Re=100 u-velocity column was cross-checked against a second,
independent web extraction of the same table -- the two agreed to the
digits both reported. The Re=3200 u-velocity column in the first gist
contains an obviously-wrong outlier (-0.86636 at y=0.4531, off by roughly
20x from its neighbors) and the second gist's own documentation flags a
separate suspect point in its Re=400 column (x=0.9063) -- both are exactly
why re-hosted transcriptions of a 1982 paper are treated as needing
independent cross-checking rather than trusted on sight. Neither flagged
value is used here: only Re=100, 400 and 1000 are compared, none of which
were flagged as suspect in either source, and Re=100 (the column actually
cross-checked twice) is the primary comparison.

Reference: Ghia, U., Ghia, K.N., Shin, C.T. (1982), "High-Re solutions for
incompressible flow using the Navier-Stokes equations and a multigrid
method", Journal of Computational Physics, 48, 387-411.

**Result, measured rather than assumed: this solver does NOT match Ghia's
table directly with its default (regularized) lid, and refining the mesh
does not close the gap -- because the two are not solving quite the same
problem, confirmed by ablation, not just argued for.** At Re=100, n=64
gives u rms/max error 0.062/0.119 against the vertical-centerline table;
n=128 gives 0.063/0.123 -- essentially unchanged despite 4x the cells,
which rules out ordinary discretization error as the explanation (that
would shrink with refinement) and points at a difference in the *posed*
boundary condition instead.

It is one: `LidDrivenCavitySolver2D`'s lid speed defaults to
`lidVelocity * sin(pi*x/L)^2`, tapered smoothly to zero at the two top
corners to remove the corner pressure/vorticity singularity a literal
discontinuous lid has -- a real, deliberate, documented design choice
(see the class's own header), not a bug. Ghia's classical cavity uses a
uniform lid instead. The regularized-cavity literature confirms this kind
of comparison typically shows only small deviations concentrated near the
corners with otherwise good agreement -- but that literature's own
regularizations usually keep the lid at full speed over most of its width
and taper only a short zone near each corner. This solver's sin^2 profile
tapers over the *entire* width, so its average speed is exactly half of a
uniform lid's (integral of sin^2 over one half-period = 1/2 -- verified
in main() below, not asserted) -- a much larger reduction in total
momentum input than a corner-only regularization, which is the
quantitative reason the deviation measured here (systemic, ~10-20% of the
field's own range, present even at the exact centerline where the local
taper value is 1) is bigger than "small deviations near the corners".

**The ablation confirms it, decisively.** `LidDrivenCavitySolver2D` gained
a `taper_lid` construction parameter (default `True`, so every existing
result and test of this class is untouched -- see the class header)
specifically to test this: at Re=100, n=64, with `taper_lid=False` (the
classical uniform lid, Ghia's own posed problem) the errors drop from
u rms/max 0.0618/0.1187 to **0.0041/0.0087**, and v rms/max from
0.0318/0.0531 to **0.0095/0.0175** -- a 3x to 15x reduction, landing under
1% of the field's own range. The gap was the lid profile, not a
discretization or physics defect, and now that is measured rather than
inferred.

**What this does and does not establish.** It does not mean the *default*
solver is wrong -- mass and momentum conservation, order-of-accuracy and
the Fourier/exact-solution comparisons elsewhere in this project's test
suite already check correctness on this same numerics, and none of that
is touched by which lid profile is chosen; the regularized lid remains
the right default for every use of this class that is not specifically
trying to reproduce Ghia's own classical problem. What it establishes is
that this engine's incompressible solver, given the same boundary
condition Ghia used, reproduces their published result to within ~1% --
the external-credibility claim this validation exercise set out to make,
now actually made rather than approximated.

Run: python python/research/ghia_1982_validation.py
"""

import sys
import time

sys.path.insert(0, "python")
import aether

# Table I: u-velocity along the vertical centerline (x=0.5), by y.
GHIA_U_Y = [1.0000, 0.9766, 0.9688, 0.9609, 0.9531, 0.8516, 0.7344, 0.6172,
            0.5000, 0.4531, 0.2813, 0.1719, 0.1016, 0.0703, 0.0625, 0.0547, 0.0000]
GHIA_U = {
    100: [1.00000, 0.84123, 0.78871, 0.73722, 0.68717, 0.23151, 0.00332, -0.13641,
          -0.20581, -0.21090, -0.15662, -0.10150, -0.06434, -0.04775, -0.04192, -0.03717, 0.00000],
    400: [1.00000, 0.75837, 0.68439, 0.61756, 0.55892, 0.29093, 0.16256, 0.02135,
          -0.11477, -0.17119, -0.32726, -0.24299, -0.14612, -0.10338, -0.09266, -0.08186, 0.00000],
    1000: [1.00000, 0.65928, 0.57492, 0.51117, 0.46604, 0.33304, 0.18719, 0.05702,
           -0.06080, -0.10648, -0.27805, -0.38289, -0.29730, -0.22220, -0.20196, -0.18109, 0.00000],
}

# Table II: v-velocity along the horizontal centerline (y=0.5), by x.
GHIA_V_X = [1.00000, 0.9688, 0.9609, 0.9531, 0.9453, 0.9063, 0.8594, 0.8047,
            0.5000, 0.2344, 0.2266, 0.1563, 0.0938, 0.0781, 0.0703, 0.0625, 0.0000]
GHIA_V = {
    100: [0.00000, -0.05906, -0.07391, -0.08864, -0.10313, -0.16914, -0.22445, -0.24533,
          0.05454, 0.17527, 0.17507, 0.16077, 0.12317, 0.10890, 0.10091, 0.09233, 0.00000],
    400: [0.00000, -0.12146, -0.15663, -0.19254, -0.22847, -0.23827, -0.44993, -0.38598,
          0.05186, 0.30174, 0.30203, 0.28124, 0.22965, 0.20920, 0.19713, 0.18360, 0.00000],
    1000: [0.00000, -0.21388, -0.27669, -0.33714, -0.39188, -0.51500, -0.42665, -0.31966,
           0.02526, 0.32235, 0.33075, 0.37095, 0.32627, 0.30353, 0.29012, 0.27485, 0.00000],
}

LENGTH = 1.0
LID_VELOCITY = 1.0


def interp(xs, ys, query):
    """Piecewise-linear interpolation over already-sorted xs (ascending)."""
    if query <= xs[0]:
        return ys[0]
    if query >= xs[-1]:
        return ys[-1]
    for k in range(1, len(xs)):
        if xs[k] >= query:
            t = (query - xs[k - 1]) / (xs[k] - xs[k - 1])
            return ys[k - 1] + t * (ys[k] - ys[k - 1])
    return ys[-1]


def converge(solver, n, tol=1e-6, max_steps=400000, check_every=200):
    dt = solver.stable_time_step()

    def snapshot():
        return [(solver.u(i, j), solver.v(i, j)) for i in range(n) for j in range(n)]

    prev = snapshot()
    steps = 0
    while steps < max_steps:
        for _ in range(check_every):
            solver.step(dt)
        steps += check_every
        cur = snapshot()
        scale = max(max(abs(u) + abs(v) for u, v in cur), 1e-12)
        change = max(abs(a[0] - b[0]) + abs(a[1] - b[1]) for a, b in zip(cur, prev)) / scale
        prev = cur
        if change <= tol:
            return steps, dt, change, True
    return steps, dt, change, False


def compare(n, reynolds, taper_lid=True):
    nu = LID_VELOCITY * LENGTH / reynolds
    scheme = aether.LidDrivenCavitySolver2D.ConvectionScheme.LIMITED_LINEAR_UPWIND
    solver = aether.LidDrivenCavitySolver2D(n, n, LENGTH, LENGTH, nu, LID_VELOCITY, scheme, taper_lid)
    t0 = time.time()
    steps, dt, change, converged = converge(solver, n)
    elapsed = time.time() - t0

    h = LENGTH / n
    # u along the vertical centerline x=0.5: nearest column of cells,
    # sampled at each cell's own y center.
    i_mid = n // 2
    xs_u = [(j + 0.5) * h for j in range(n)]
    us = [solver.u(i_mid, j) for j in range(n)]
    # The y=0/y=1 rows are excluded from scoring: they are Ghia's *boundary
    # condition* (u=0 at the floor, u=1 at the lid, by construction, not a
    # predicted quantity), while this solver is collocated and has no value
    # defined exactly at the wall -- its nearest cell center sits half a
    # cell inward. Comparing that cell center against the wall's own
    # prescribed value would score a structural difference between
    # collocated and boundary-defined grids, not a physical one.
    u_pairs = [(y, interp(xs_u, us, y), g) for y, g in zip(GHIA_U_Y, GHIA_U[reynolds])
               if 0.0 < y < 1.0]
    u_errors = [abs(mine - ghia) for _, mine, ghia in u_pairs]

    j_mid = n // 2
    xs_v = [(i + 0.5) * h for i in range(n)]
    vs = [solver.v(i, j_mid) for i in range(n)]
    v_pairs = [(x, interp(xs_v, vs, x), g) for x, g in zip(GHIA_V_X, GHIA_V[reynolds])
               if 0.0 < x < 1.0]
    v_errors = [abs(mine - ghia) for _, mine, ghia in v_pairs]

    def rms(values):
        return (sum(v * v for v in values) / len(values)) ** 0.5

    return {
        "n": n, "reynolds": reynolds, "steps": steps, "elapsed": elapsed,
        "converged": converged, "change": change,
        "u_rms": rms(u_errors), "u_max": max(u_errors),
        "v_rms": rms(v_errors), "v_max": max(v_errors),
    }


def main():
    print("LidDrivenCavitySolver2D vs Ghia, Ghia & Shin (1982)")
    print(f"{'Re':>6} {'n':>4} {'lid':>8} {'passos':>8} {'tempo':>7} {'convergiu':>10} "
          f"{'u rms':>10} {'u max':>10} {'v rms':>10} {'v max':>10}")
    for reynolds, n in [(100, 64), (100, 128), (400, 64), (1000, 64)]:
        result = compare(n, reynolds)
        print(f"{result['reynolds']:>6} {result['n']:>4} {'tapered':>8} {result['steps']:>8} "
              f"{result['elapsed']:>6.1f}s {str(result['converged']):>10} "
              f"{result['u_rms']:>10.4f} {result['u_max']:>10.4f} "
              f"{result['v_rms']:>10.4f} {result['v_max']:>10.4f}")

    # The ablation: same solver, same mesh, only the lid's own boundary
    # condition changed to the classical uniform one -- LidDrivenCavitySolver2D's
    # taper_lid=False, added specifically for this comparison. If the taper
    # is indeed what separates this solver's result from Ghia's table, the
    # errors below should drop sharply relative to the "tapered" row at the
    # same (Re, n) above; if they do not, the taper is not the (whole)
    # explanation and the module docstring's open question stands.
    # n=64 only: the tapered comparison above already showed n=64 vs n=128
    # give essentially the same error, so a second, much more expensive
    # n=128 point would add runtime without adding evidence either way.
    print("\nablacao: mesmo solver, lid uniforme (taper_lid=False) em vez do lid regularizado")
    for reynolds, n in [(100, 64)]:
        result = compare(n, reynolds, taper_lid=False)
        print(f"{result['reynolds']:>6} {result['n']:>4} {'uniform':>8} {result['steps']:>8} "
              f"{result['elapsed']:>6.1f}s {str(result['converged']):>10} "
              f"{result['u_rms']:>10.4f} {result['u_max']:>10.4f} "
              f"{result['v_rms']:>10.4f} {result['v_max']:>10.4f}")

    # See this file's own module docstring for the finding this number
    # explains: LidDrivenCavitySolver2D's lid speed is lidVelocity*sin(pi*x)^2,
    # not a uniform lidVelocity -- its average over the lid's own width is
    # this integral, a closed-form fact (integral of sin^2 over one full
    # period-half is exactly half the domain length), not a simulation
    # result. Printed here so the claim in the docstring is checked by the
    # same run that produced the numbers it explains, not left as an
    # unverified aside.
    import math
    samples = 100000
    mean_taper = sum(math.sin(math.pi * (i + 0.5) / samples) ** 2 for i in range(samples)) / samples
    print(f"\nfator de taper medio do lid (sin^2(pi*x), x em [0,1]): {mean_taper:.4f} "
          "(exato: 0.5) -- metade do momento de um lid uniforme")


if __name__ == "__main__":
    main()
