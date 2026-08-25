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


def _background_lattice_points(box_min: Vector3, box_max: Vector3, spacing: float,
                                exclude_center: Vector3, exclude_radius: float,
                                seed: int = 17, jitter_fraction: float = 0.25) -> list[Vector3]:
    """A lattice of background points spanning the box, at roughly `spacing`
    apart, with the object's own neighbourhood carved out (see
    `exclude_radius`'s caller).

    **Reuses, rather than invents, the box-lattice convention this project
    already validated everywhere else**: the six boundary layers (i, j or k
    at 0 or n) sit *exactly* on the box's planes -- unperturbed, because
    `classify_boundary_face`'s plane-membership test depends on that -- while
    every interior point is jittered by up to `jitter_fraction` of its local
    cell spacing, the same jitter `build_jittered_lattice()`
    (python/tests/test_unstructured_bindings.py, used by every unstructured
    solver test since DIVIDA_TECNICA.md's whole 4.x investigation) applies
    for the same reason: an unperturbed regular lattice tetrahedralizes into
    exact co-planar/co-spherical ties, which is Delaunay's known degenerate
    case, not a representative test of a real mesh.

    `exclude_radius` is a sphere around `exclude_center`, not a distance to
    the object's actual surface -- cheap, and correct for a roughly round or
    star-shaped object, but conservative-to-wasteful for a very elongated
    one (a long thin object's "waist" sits far from its own bounding
    sphere's edge, so background points get excluded there too, coarser
    than necessary). Fixing that needs distance-to-surface, not
    distance-to-bounding-sphere -- a real next step, not attempted here.
    """
    rng = random.Random(seed)
    nx = max(2, round((box_max.x - box_min.x) / spacing))
    ny = max(2, round((box_max.y - box_min.y) / spacing))
    nz = max(2, round((box_max.z - box_min.z) / spacing))
    points: list[Vector3] = []
    for i in range(nx + 1):
        x = box_min.x + (box_max.x - box_min.x) * i / nx
        if 0 < i < nx:
            x += rng.uniform(-jitter_fraction, jitter_fraction) * (box_max.x - box_min.x) / nx
        for j in range(ny + 1):
            y = box_min.y + (box_max.y - box_min.y) * j / ny
            if 0 < j < ny:
                y += rng.uniform(-jitter_fraction, jitter_fraction) * (box_max.y - box_min.y) / ny
            for k in range(nz + 1):
                z = box_min.z + (box_max.z - box_min.z) * k / nz
                if 0 < k < nz:
                    z += rng.uniform(-jitter_fraction, jitter_fraction) * (box_max.z - box_min.z) / nz
                dx, dy, dz = x - exclude_center.x, y - exclude_center.y, z - exclude_center.z
                if math.sqrt(dx * dx + dy * dy + dz * dz) < exclude_radius:
                    continue
                points.append(Vector3(x, y, z))
    return points


def mesh_flow_around_object(triangle_mesh: TriangleMesh, margin: float = 3.0,
                             background_spacing: float | None = None,
                             background_coarsening: float = 2.0,
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

    **Background points, not just box corners.** The first version of this
    function only added the box's 8 corners besides the object's own
    surface -- correct topologically, but it meant almost the entire empty
    bulk of the box became a handful of enormous tetrahedra: measured on a
    real test case, an 80x ratio between the smallest and largest cell,
    which was enough on its own to blow the solver up (divergence 3.9e6 in
    50 steps) with mesh-quality diagnostics (non-orthogonality, deficient
    stencils) that looked completely healthy -- this was never the
    DIVIDA_TECNICA.md 4.3 instability, just an under-resolved bulk. A first
    fix attempt (refining after the fact via `insert_steiner_point` on the
    largest cells) corrupted 5 of 20 already-recovered object facets, because
    Steiner insertion has no notion of a protected wall. `background_spacing`
    (default: `background_coarsening` times the object surface's own average
    triangle edge length) now seeds a jittered background lattice *before*
    `tetrahedralize()` runs at all, via `_background_lattice_points` --
    avoiding the problem instead of repairing it after the fact. Points
    within `exclude_margin` times the object's bounding-sphere radius of its
    centroid are dropped, so the background lattice never crowds the
    object's own (finer) surface resolution.

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
