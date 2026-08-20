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


def mesh_flow_around_object(triangle_mesh: TriangleMesh, margin: float = 3.0) -> FlowDomain:
    """Builds a fluid-domain `TetrahedralMesh` for flow around `triangle_mesh`
    inside an axis-aligned bounding box expanded by `margin` times the
    object's own bounding-box size on every side.

    `triangle_mesh` must already be watertight (call `weld_vertices()` first
    if it came straight from `load_stl`/`load_obj` and hasn't been welded) --
    checked explicitly rather than left to fail confusingly deep inside
    `remove_region`. The object's own vertex indices are reused unchanged as
    tetrahedralization point indices (points 0..N-1), so its triangles can be
    handed to `recover_facets`/`remove_region` without any index remapping.

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

    for x in (box_min.x, box_max.x):
        for y in (box_min.y, box_max.y):
            for z in (box_min.z, box_max.z):
                tetra.add_point(x, y, z)

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
