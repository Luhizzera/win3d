"""Geometry -> volume mesh -> solver glue for external flow around an
imported object (Roadmap item A1: "algo palpável e utilizável em escopo
real", fase A).

**What this closes.** Before this module, three pieces existed in complete
isolation: STL/OBJ import (`aether.load_stl`/`load_obj`), boundary-recovery
tetrahedralization (`DelaunayTetrahedralization3D.recover_facets`/
`remove_region`), and the unstructured Navier-Stokes solver
(`UnstructuredCavitySolver3D`). Nothing chained them -- confirmed by grep
before writing this: `load_stl`/`load_obj` were used nowhere outside their
own bindings, and the only script exercising `UnstructuredCavitySolver3D`
built its mesh exclusively from a synthetic jittered lattice, never an
imported surface.

**Scope chosen deliberately: external flow, not internal (duct) flow.**
`remove_region()`'s own shape -- flood-fill from an interior seed point,
stopping at recovered walls -- maps directly onto "carve the object's solid
interior out of a bounding box", i.e. external aerodynamics (flow around a
car, a building, an airfoil). It does *not* map onto internal duct/pipe
flow (the interior of a possibly non-convex closed surface, discarding
tetrahedra *outside* the recovered walls instead of inside a seed region) --
that is a real, different capability, not attempted here.

**Scope chosen deliberately: closed domain, not inlet/outlet.** The mesh
this module builds is a plain box with the object as an interior wall --
callers are expected to build a *closed* `UnstructuredCavitySolver3D`
(no `is_outlet`), driving flow via a moving wall (like a lid) on one box
face, the way `driving_wall_velocity()` below does. DIVIDA_TECNICA.md 4.3's
coupling-correction fix is proven only for closed domains; a domain with an
outlet still has the open GMRES-stagnation gap documented there. Extending
this pipeline to inlet/outlet domains is future work gated on that fix, not
on anything in this module.
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass

from aether_core_py import Vector3
from aether_geometry_py import TriangleMesh
from aether_mesh_py import DelaunayTetrahedralization3D, TetrahedralMesh


class MeshGenerationError(RuntimeError):
    """Raised when the geometry->volume-mesh chain cannot produce a usable
    fluid-domain mesh, instead of silently returning something wrong.

    Every case this covers was a specific, checkable failure mode of
    `remove_region`/`recover_facets`, not a generic catch-all: an
    unrecoverable facet, a seed point that landed outside every tetrahedron
    (so nothing was removed), or a carved volume that does not plausibly
    match "box minus object".
    """


@dataclass(frozen=True)
class FlowDomain:
    """The result of `mesh_flow_around_object`: a fluid-domain
    `TetrahedralMesh` plus enough geometric bookkeeping for a caller to
    classify its own boundary faces (see `classify_boundary_face`) without
    needing to touch the tetrahedralizer again.
    """

    mesh: TetrahedralMesh
    box_min: Vector3
    box_max: Vector3
    unrecovered_facet_count: int
    carved_volume: float
    expected_volume: float
    background_point_count: int


def _bounding_box(triangle_mesh: TriangleMesh) -> tuple[Vector3, Vector3]:
    lo = triangle_mesh.vertex(0)
    hi = triangle_mesh.vertex(0)
    lo = Vector3(lo.x, lo.y, lo.z)
    hi = Vector3(hi.x, hi.y, hi.z)
    for i in range(1, triangle_mesh.vertex_count()):
        v = triangle_mesh.vertex(i)
        lo = Vector3(min(lo.x, v.x), min(lo.y, v.y), min(lo.z, v.z))
        hi = Vector3(max(hi.x, v.x), max(hi.y, v.y), max(hi.z, v.z))
    return lo, hi


def _average_surface_edge_length(triangle_mesh: TriangleMesh) -> float:
    """A cheap proxy for the object surface's own resolution: sqrt(area /
    triangle_count) is the edge length an equilateral triangle of the mesh's
    average triangle area would have. Exact edge lengths would need walking
    every triangle's three edges; this is the same order of magnitude for
    any reasonably-shaped mesh and costs one call each of two quantities
    already computed by TriangleMesh itself.
    """
    n = triangle_mesh.triangle_count()
    if n == 0:
        return 1.0
    return math.sqrt(triangle_mesh.surface_area() / n)


def _background_lattice_points(box_min: Vector3, box_max: Vector3, finest_spacing: float,
                                exclude_center: Vector3, exclude_radius: float,
                                max_level: int = 0, seed: int = 17,
                                jitter_fraction: float = 0.25) -> list[Vector3]:
    """A **graded** lattice of background points spanning the box: spacing
    `finest_spacing` in the shell just outside the object, doubling with
    every doubling of distance from it, up to `2**max_level` times coarser
    in the far field. The object's own neighbourhood (inside
    `exclude_radius`) is left empty for its surface points to own.

    **`max_level` defaults to 0, which means no grading at all -- and that
    default was measured, not assumed.** Grading was built to attack an
    apparent seam where the object's surface resolution met a background
    lattice of a different density (worst-to-best cell volume ratio ~97x,
    max non-orthogonality 2.76). A direct A/B on the same geometry and
    spacing settled it: graded (max_level=3) gives 886 cells from 168
    background points with volume ratio 965 and non-orthogonality 2.52,
    while uniform gives 3493 cells from 722 points with ratio **88** and
    non-orthogonality **2.06**. Grading is cheaper but measurably *worse*
    on both quality metrics, because each decimation threshold introduces
    its own abrupt density jump -- trading one seam for `max_level` of
    them. Both meshes turned out to run the solver correctly once an
    unrelated bug in the driving boundary condition was fixed (see
    `driving_wall_velocity`), so the seam was never the blocker it looked
    like; uniform is the default purely because better quality metrics
    leave more margin. Grading is kept, and kept documented, because
    coarsening the far field is the right *idea* -- it just needs smooth
    size variation (Delaunay refinement against a size function) rather
    than lattice decimation to pay off.

    **How the grading stays consistent.** Rather than generating each shell
    independently (which would put unrelated point lattices next to each
    other and create a new seam at every shell boundary), one finest-level
    lattice is defined over the whole box and then *decimated*: a point at
    index (i, j, k) survives only if all three indices are multiples of
    2**level for its own distance-derived level. Every coarse point is
    therefore also a fine-lattice point, so consecutive levels stay
    aligned -- the same nesting an octree gives. `nx/ny/nz` are rounded up
    to a multiple of `2**max_level` so that the far index (i == nx) is a
    multiple of every level's stride and the box's own boundary planes
    survive decimation.

    **Reuses, rather than invents, the box-lattice convention this project
    already validated everywhere else**: the six boundary layers (i, j or k
    at 0 or n) sit *exactly* on the box's planes -- unperturbed, because
    `classify_boundary_face`'s plane-membership test depends on that --
    while every interior point is jittered, here by `jitter_fraction` of
    its *own level's* spacing rather than the finest one, so coarse regions
    get proportionally the same perturbation. The jitter itself is the same
    idea `build_jittered_lattice()`
    (python/tests/test_unstructured_bindings.py, used by every unstructured
    solver test since DIVIDA_TECNICA.md's whole 4.x investigation) applies
    for the same reason: an unperturbed regular lattice tetrahedralizes
    into exact co-planar/co-spherical ties, which is Delaunay's known
    degenerate case, not a representative test of a real mesh.

    `exclude_radius` is a sphere around `exclude_center`, not a distance to
    the object's actual surface -- cheap, and correct for a roughly round or
    star-shaped object, but conservative-to-wasteful for a very elongated
    one (a long thin object's "waist" sits far from its own bounding
    sphere's edge, so background points get excluded there too, coarser
    than necessary). Fixing that needs distance-to-surface, not
    distance-to-bounding-sphere -- a real next step, not attempted here.
    """
    rng = random.Random(seed)
    stride = 2 ** max_level

    def _axis_count(extent: float) -> int:
        n = max(1, math.ceil(extent / finest_spacing))
        return int(math.ceil(n / stride) * stride)

    ex, ey, ez = box_max.x - box_min.x, box_max.y - box_min.y, box_max.z - box_min.z
    nx, ny, nz = _axis_count(ex), _axis_count(ey), _axis_count(ez)
    hx, hy, hz = ex / nx, ey / ny, ez / nz

    points: list[Vector3] = []
    for i in range(nx + 1):
        for j in range(ny + 1):
            for k in range(nz + 1):
                x0 = box_min.x + ex * i / nx
                y0 = box_min.y + ey * j / ny
                z0 = box_min.z + ez * k / nz
                dx, dy, dz = x0 - exclude_center.x, y0 - exclude_center.y, z0 - exclude_center.z
                radius = math.sqrt(dx * dx + dy * dy + dz * dz)
                if radius < exclude_radius:
                    continue
                # Level 0 hugs the exclusion sphere and each level out
                # covers a doubling of radius; clamped so the far field
                # never gets coarser than max_level.
                level = int(math.floor(math.log2(radius / exclude_radius))) if radius > exclude_radius else 0
                level = max(0, min(max_level, level))
                step = 2 ** level
                if i % step or j % step or k % step:
                    continue
                x, y, z = x0, y0, z0
                if 0 < i < nx:
                    x += rng.uniform(-jitter_fraction, jitter_fraction) * hx * step
                if 0 < j < ny:
                    y += rng.uniform(-jitter_fraction, jitter_fraction) * hy * step
                if 0 < k < nz:
                    z += rng.uniform(-jitter_fraction, jitter_fraction) * hz * step
                points.append(Vector3(x, y, z))
    return points


def mesh_flow_around_object(triangle_mesh: TriangleMesh, margin: float = 3.0,
                             background_spacing: float | None = None,
                             background_coarsening: float = 1.0,
                             exclude_margin: float = 1.15) -> FlowDomain:
    """Builds a fluid-domain `TetrahedralMesh` for flow around `triangle_mesh`
    inside an axis-aligned bounding box expanded by `margin` times the
    object's own bounding-box size on every side.

    `triangle_mesh` must already be watertight (call `weld_vertices()` first
    if it came straight from `load_stl`/`load_obj` and hasn't been welded) --
    checked explicitly rather than left to fail confusingly deep inside
    `remove_region`. The object's own vertex indices are reused unchanged as
    tetrahedralization point indices (points 0..N-1), so its triangles can be
    handed to `recover_facets`/`remove_region` without any index remapping.

    **Background points, not just box corners, and graded.** The first
    version of this function only added the box's 8 corners besides the
    object's own surface -- correct topologically, but it meant almost the
    entire empty bulk of the box became a handful of enormous tetrahedra:
    measured on a real test case, an 80x ratio between the smallest and
    largest cell, which was enough on its own to blow the solver up
    (divergence 3.9e6 in 50 steps) with mesh-quality diagnostics
    (non-orthogonality, deficient stencils) that looked completely healthy
    -- this was never the DIVIDA_TECNICA.md 4.3 instability, just an
    under-resolved bulk. A first fix attempt (refining after the fact via
    `insert_steiner_point` on the largest cells) corrupted 5 of 20
    already-recovered object facets, because Steiner insertion has no
    notion of a protected wall. Seeding the background *before*
    `tetrahedralize()` avoids that problem entirely rather than repairing
    it afterwards, and `_background_lattice_points` grades the result so
    the background meets the object's own surface resolution instead of
    colliding with it -- see that function for the measurements behind
    both decisions.

    `background_spacing` (default: `background_coarsening` times the object
    surface's own average triangle edge length) is the *finest* spacing,
    used in the shell immediately outside the object; the far field is
    automatically coarser. Points within `exclude_margin` times the
    object's bounding-sphere radius of its centroid are left out, so the
    object's own surface points own that region.

    Raises `MeshGenerationError` if any facet could not be recovered, if the
    interior seed point (the object's own vertex centroid -- correct for a
    convex or star-shaped object, the same assumption
    `TriangleMesh.reorient_normals_outward()` already documents elsewhere in
    this engine) did not land inside any tetrahedron, or if the carved volume
    is not within 1% of (box volume - object volume) -- the cheapest
    catchable signal that recovery silently went wrong somewhere else in the
    domain, since `remove_region`'s flood-fill has no way to know on its own
    whether it stopped at the intended walls.
    """
    if not triangle_mesh.is_watertight():
        raise MeshGenerationError(
            "mesh_flow_around_object: the object mesh is not watertight -- "
            "call weld_vertices() first, or the surface genuinely has a hole."
        )

    box_lo, box_hi = _bounding_box(triangle_mesh)
    size = Vector3(box_hi.x - box_lo.x, box_hi.y - box_lo.y, box_hi.z - box_lo.z)
    box_min = Vector3(box_lo.x - margin * size.x, box_lo.y - margin * size.y,
                       box_lo.z - margin * size.z)
    box_max = Vector3(box_hi.x + margin * size.x, box_hi.y + margin * size.y,
                       box_hi.z + margin * size.z)

    tetra = DelaunayTetrahedralization3D()
    centroid = Vector3(0.0, 0.0, 0.0)
    for i in range(triangle_mesh.vertex_count()):
        v = triangle_mesh.vertex(i)
        tetra.add_point(v.x, v.y, v.z)
        centroid = Vector3(centroid.x + v.x, centroid.y + v.y, centroid.z + v.z)
    centroid = Vector3(centroid.x / triangle_mesh.vertex_count(),
                        centroid.y / triangle_mesh.vertex_count(),
                        centroid.z / triangle_mesh.vertex_count())

    bbox_center = Vector3((box_lo.x + box_hi.x) / 2.0, (box_lo.y + box_hi.y) / 2.0,
                           (box_lo.z + box_hi.z) / 2.0)
    bbox_half_diagonal = 0.5 * math.sqrt(size.x ** 2 + size.y ** 2 + size.z ** 2)
    exclude_radius = exclude_margin * bbox_half_diagonal

    if background_spacing is None:
        background_spacing = background_coarsening * _average_surface_edge_length(triangle_mesh)

    background_points = _background_lattice_points(box_min, box_max, background_spacing,
                                                     bbox_center, exclude_radius)
    for p in background_points:
        tetra.add_point(p.x, p.y, p.z)

    tetra.tetrahedralize()

    facets = [triangle_mesh.triangle(i) for i in range(triangle_mesh.triangle_count())]
    recovery = tetra.recover_facets(facets)
    if recovery.unrecovered:
        raise MeshGenerationError(
            f"mesh_flow_around_object: {len(recovery.unrecovered)} of "
            f"{len(facets)} object facets could not be recovered (Schonhardt-"
            "style obstruction or badly-conditioned input geometry) -- see "
            "DelaunayTetrahedralization3D.recover_facets' own doc comment. "
            "Refine or clean the input mesh before retrying."
        )

    removed = tetra.remove_region(centroid, recovery.recovered_facets)
    if removed == 0:
        raise MeshGenerationError(
            "mesh_flow_around_object: the object's vertex centroid did not "
            "land inside any tetrahedron -- it is likely not star-shaped "
            "(reorient_normals_outward()'s own documented assumption), so "
            "the centroid sits outside the solid. Supply a known-interior "
            "seed point directly for a non-star-shaped object."
        )

    mesh = TetrahedralMesh.from_tetrahedralization(tetra)

    box_volume = (box_max.x - box_min.x) * (box_max.y - box_min.y) * (box_max.z - box_min.z)
    expected_volume = box_volume - triangle_mesh.volume()
    carved_volume = mesh.total_volume()
    if abs(carved_volume - expected_volume) > 0.01 * expected_volume:
        raise MeshGenerationError(
            f"mesh_flow_around_object: carved fluid volume {carved_volume:.6g} "
            f"does not match box-minus-object {expected_volume:.6g} within 1% "
            "-- remove_region likely leaked past a wall or removed the wrong "
            "region. Not returning a mesh that fails its own volume check."
        )

    return FlowDomain(
        mesh=mesh,
        box_min=box_min,
        box_max=box_max,
        unrecovered_facet_count=len(recovery.unrecovered),
        carved_volume=carved_volume,
        expected_volume=expected_volume,
        background_point_count=len(background_points),
    )


def classify_boundary_face(position: Vector3, domain: FlowDomain, tol: float = 1e-6) -> str:
    """Which boundary this face belongs to: one of "x_min", "x_max", "y_min",
    "y_max", "z_min", "z_max" (a face of the surrounding box) or "object"
    (a face of the carved-out obstacle -- anywhere not on the box's own
    planes, by construction, since the object sits strictly inside the box
    by `margin`).

    Written to plug directly into `UnstructuredCavitySolver3D`'s
    `wall_velocity(position)` selector: a caller matches on this string to
    decide which box face drives the flow and confirms every other face
    (including "object") gets a no-slip zero-velocity condition.
    """
    if abs(position.x - domain.box_min.x) < tol:
        return "x_min"
    if abs(position.x - domain.box_max.x) < tol:
        return "x_max"
    if abs(position.y - domain.box_min.y) < tol:
        return "y_min"
    if abs(position.y - domain.box_max.y) < tol:
        return "y_max"
    if abs(position.z - domain.box_min.z) < tol:
        return "z_min"
    if abs(position.z - domain.box_max.z) < tol:
        return "z_max"
    return "object"


# Outward axis normal of each box face, used to enforce tangency below.
_FACE_NORMALS = {
    "x_min": (-1.0, 0.0, 0.0), "x_max": (1.0, 0.0, 0.0),
    "y_min": (0.0, -1.0, 0.0), "y_max": (0.0, 1.0, 0.0),
    "z_min": (0.0, 0.0, -1.0), "z_max": (0.0, 0.0, 1.0),
}


def driving_wall_velocity(domain: FlowDomain, face: str = "z_max",
                           direction: tuple[float, float, float] = (1.0, 0.0, 0.0),
                           speed: float = 1.0):
    """Builds the `wall_velocity(position)` callable
    `UnstructuredCavitySolver3D` expects: one box face slides in its own
    plane at `speed`, tapered smoothly to zero at that face's four edges,
    and every other boundary -- the five remaining box faces and the
    object itself -- is a stationary no-slip wall.

    **This exists because getting it wrong is easy and the failure is
    dramatic.** A closed domain has no outlet, so the prescribed wall
    motion must inject no mass: the driving velocity has to be *tangential*
    to the face carrying it. Prescribing a velocity with a component along
    that face's own normal instead blows fluid into a sealed box, which is
    not a hard case for the solver so much as an unsolvable one -- mass
    accumulates with nowhere to go. Measured while building this module,
    driving the `x_min` face (normal x) in the x direction: peak velocity
    reached 119x the driving speed and face divergence 1e2, on a mesh whose
    spectral radius around rest was a perfectly stable 0.245 and whose rest
    state held at exactly zero. The mesh was never the problem. The solver
    reported the real cause the whole time through `net_boundary_flux()`,
    which sat at -26 instead of 0; with tangency enforced the same mesh
    gives net flux exactly 0.00e+00, divergence ~1e-8, and the expected
    primary-vortex topology.

    Raises ValueError rather than silently projecting out the offending
    component: a caller who asked for a normal-direction push has a
    misconception worth surfacing, not a typo worth quietly repairing.

    The taper is the same `sin^2` regularization
    `LidDrivenCavitySolver2D` uses, for the same documented reason: a lid
    that slides at full speed right up to where it meets a stationary wall
    is discontinuous there, and that discontinuity is a genuine pressure
    singularity of the continuous problem, not a discretization artifact.
    """
    if face not in _FACE_NORMALS:
        raise ValueError(f"driving_wall_velocity: unknown face {face!r}; "
                          f"expected one of {sorted(_FACE_NORMALS)}")
    nx, ny, nz = _FACE_NORMALS[face]
    dx, dy, dz = direction
    along_normal = dx * nx + dy * ny + dz * nz
    magnitude = math.sqrt(dx * dx + dy * dy + dz * dz)
    if magnitude == 0.0:
        raise ValueError("driving_wall_velocity: direction must be nonzero")
    if abs(along_normal) > 1e-12 * magnitude:
        raise ValueError(
            f"driving_wall_velocity: direction {direction} is not tangential to face "
            f"{face!r} (normal {(nx, ny, nz)}). A closed domain cannot absorb the mass "
            "a normal-direction wall velocity injects -- see this function's own doc "
            "comment for what that failure looks like."
        )
    ux, uy, uz = dx / magnitude * speed, dy / magnitude * speed, dz / magnitude * speed

    # The two in-plane axes of this face, and the box extent along each,
    # for the taper's normalized coordinates.
    axes = [a for a in range(3) if (nx, ny, nz)[a] == 0.0]
    lo = (domain.box_min.x, domain.box_min.y, domain.box_min.z)
    hi = (domain.box_max.x, domain.box_max.y, domain.box_max.z)

    def wall_velocity(position: Vector3) -> Vector3:
        if classify_boundary_face(position, domain) != face:
            return Vector3(0.0, 0.0, 0.0)
        coords = (position.x, position.y, position.z)
        taper = 1.0
        for a in axes:
            span = hi[a] - lo[a]
            t = (coords[a] - lo[a]) / span if span != 0.0 else 0.5
            taper *= math.sin(math.pi * min(max(t, 0.0), 1.0)) ** 2
        return Vector3(ux * taper, uy * taper, uz * taper)

    return wall_velocity
