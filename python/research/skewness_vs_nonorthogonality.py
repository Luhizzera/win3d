"""DIVIDA_TECNICA.md 3.3: separate skewness from non-orthogonality.

The item's own measurement used geometrically graded lattices, where the
two rise *together* (non-orthogonality 1.63 -> 4.90 -> 15.08 alongside
skewness 1.15 -> 3.17 -> 10.05), so it could not say which one costs the
order. This builds two families that move one while holding the other
roughly fixed:

  SHEAR   -- an affine shear maps every cell to a congruent parallelepiped.
             The centroid-to-centroid line still crosses each face at its
             centre (skewness stays low) while the face normal tilts away
             from that line (non-orthogonality rises).

  GRADING -- a stretched-but-axis-aligned lattice keeps face normals
             axis-aligned (non-orthogonality stays low) while cells of
             different sizes put the crossing point off-centre (skewness
             rises).

Both are then measured the same way item 3.2 established: a manufactured
solution, whose observed order is the number under test.

**Result, and its ceiling.** The shear family isolates cleanly up to
strength 1.2: skewness stays around 0.4-1.0 while non-orthogonality climbs
1.54 -> 9.43, and the observed order stays at 2.06-2.15. Past that the
tetrahedralization stops preserving the affine structure -- at strength
1.8-2.4 skewness climbs to 2.2-2.4 as well -- so the two stop being
separable and the order falls with them. The grading family never
separated them at all (non-orthogonality rose to 16.33), so it measures
nothing this script set out to measure and is kept only to show that.

Run from the repository root:  python python/research/skewness_vs_nonorthogonality.py
"""
import math
import sys

import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."  , "python"))
import aether


def lattice(n, mode, strength, seed=17):
    """A unit-cube point lattice, distorted one of two ways.

    Boundary points stay on their planes in both cases so the domain is
    still the unit cube and the boundary selectors below stay valid.
    """
    import random
    rng = random.Random(seed)
    tets = aether.DelaunayTetrahedralization3D()
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(n + 1):
                u, v, w = i / n, j / n, k / n
                if mode == "shear":
                    # x is displaced by a multiple of z: an affine shear.
                    # Interior jitter is kept small and identical across
                    # strengths so it is not what varies.
                    x = u + strength * w
                    y, z = v, w
                    # The shear pushes x outside [0,1]; renormalise so the
                    # domain stays a box (the selectors below use planes).
                    x = (x - 0.0) / (1.0 + strength)
                elif mode == "grade":
                    # Geometric stretching along x only: axis-aligned, so
                    # face normals stay axis-aligned.
                    x = (math.exp(strength * u) - 1.0) / (math.exp(strength) - 1.0) if strength > 0 else u
                    y, z = v, w
                else:
                    x, y, z = u, v, w
                if 0 < i < n and 0 < j < n and 0 < k < n:
                    jitter = 0.12 / n
                    x += rng.uniform(-jitter, jitter)
                    y += rng.uniform(-jitter, jitter)
                    z += rng.uniform(-jitter, jitter)
                tets.add_point(x, y, z)
    tets.tetrahedralize()
    return aether.TetrahedralMesh.from_tetrahedralization(tets)


def skewness(mesh):
    """Peak face skewness: the distance from a face's centroid to where the
    owner-neighbour line crosses its plane, normalised by the face's own
    size (sqrt of its area). This is the standard definition and it is
    computed here rather than read off the engine, which does not expose
    it.
    """
    worst = 0.0
    for f in range(mesh.face_count()):
        if mesh.is_boundary_face(f):
            continue
        face = mesh.face(f)
        c0 = mesh.cell_centroid(face.owner)
        c1 = mesh.cell_centroid(face.neighbour)
        d = aether.Vector3(c1.x - c0.x, c1.y - c0.y, c1.z - c0.z)
        a = face.area_vector
        denom = a.x * d.x + a.y * d.y + a.z * d.z
        if abs(denom) < 1e-300:
            continue
        # Parameter where the segment crosses the face plane.
        fc = face.centroid
        num = a.x * (fc.x - c0.x) + a.y * (fc.y - c0.y) + a.z * (fc.z - c0.z)
        t = num / denom
        px, py, pz = c0.x + t * d.x, c0.y + t * d.y, c0.z + t * d.z
        offset = math.sqrt((px - fc.x) ** 2 + (py - fc.y) ** 2 + (pz - fc.z) ** 2)
        area = math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
        if area > 0.0:
            worst = max(worst, offset / math.sqrt(area))
    return worst


def manufactured_error(mesh):
    """phi = sin(pi x) sin(pi y) sin(pi z), whose Laplacian is
    -3 pi^2 phi -- so the source term is known exactly and the only error
    left is the discretization's own."""
    kPi = math.pi

    def exact(p):
        return math.sin(kPi * p.x) * math.sin(kPi * p.y) * math.sin(kPi * p.z)

    solver = aether.UnstructuredDiffusionSolver(mesh)
    solver.set_dirichlet_boundary(lambda p: True, exact)
    solver.set_source_term(lambda p: 3.0 * kPi * kPi * exact(p))
    sweeps = solver.solve_conjugate_gradient(20000, 1e-10, 2000)

    total = 0.0
    volume = 0.0
    for c in range(mesh.cell_count()):
        centroid = mesh.cell_centroid(c)
        error = solver.value(c) - exact(centroid)
        w = mesh.cell_volume(c)
        total += w * error * error
        volume += w
    return math.sqrt(total / volume), sweeps, solver.max_non_orthogonality()


print(f"{'familia':<10}{'forca':<8}{'n':<4}{'cel':<7}{'skew':<9}{'naoOrt':<9}{'rms':<12}{'ordem':<8}{'varr':<7}", flush=True)
for mode, strengths in (("shear", (0.0, 0.6, 1.2, 1.8, 2.4)), ("grade", (0.0, 1.5, 3.0))):
    for strength in strengths:
        previous = None
        for n in (4, 6):
            mesh = lattice(n, mode, strength)
            sk = skewness(mesh)
            rms, sweeps, nonorth = manufactured_error(mesh)
            order = ""
            if previous is not None:
                rms0, n0 = previous
                if rms > 0 and rms0 > 0:
                    order = f"{math.log(rms0 / rms) / math.log(n / n0):.2f}"
            previous = (rms, n)
            print(f"{mode:<10}{strength:<8.2f}{n:<4}{mesh.cell_count():<7}{sk:<9.2f}"
                  f"{nonorth:<9.2f}{rms:<12.4e}{order:<8}{sweeps:<7}", flush=True)
