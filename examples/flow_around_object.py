"""External flow around an imported object, end to end.

This is the whole pipeline in one file: build (or load) a closed surface,
generate a fluid-domain volume mesh around it, check that the case is
solvable *before* spending time on it, march to steady state, and write a
result ParaView can open.

Run it from the repository root after `python build.py`:

    python examples/flow_around_object.py

    python examples/flow_around_object.py meu_objeto.stl

With no argument it builds an icosahedron in memory, so the example runs
on a clean clone with no data files. With a path it loads that STL or OBJ
instead -- which is the case this example exists to demonstrate, since
everything downstream is identical either way.
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))

import aether


def build_icosahedron() -> "aether.TriangleMesh":
    """A unit-circumradius icosahedron: 12 vertices, 20 triangles, convex
    and watertight, so the example needs no data file to run."""
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    raw = [(-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
           (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
           (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1)]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
             (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
             (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    norm = math.sqrt(1.0 + phi * phi)
    mesh = aether.TriangleMesh()
    for x, y, z in raw:
        mesh.add_vertex(aether.Vector3(x / norm, y / norm, z / norm))
    for a, b, c in faces:
        mesh.add_triangle(a, b, c)
    return mesh


def load_surface(path: str) -> "aether.TriangleMesh":
    """Loads an STL or OBJ and gets it into the state the mesher needs.

    `load_stl` already welds -- STL stores three independent vertices per
    triangle with no shared connectivity, so without welding no two
    triangles share an edge and the surface looks like it is all holes.
    `load_obj` does not weld and does not need to: OBJ already stores a
    shared, indexed vertex list.
    """
    extension = os.path.splitext(path)[1].lower()
    if extension == ".obj":
        surface = aether.load_obj(path)
    else:
        surface = aether.load_stl(path)
    surface.remove_degenerate_triangles()
    return surface


def main() -> int:
    if len(sys.argv) > 1:
        surface = load_surface(sys.argv[1])
        print(f"geometria carregada de {sys.argv[1]}")
    else:
        surface = build_icosahedron()
        print("geometria: icosaedro construido em memoria (passe um .stl/.obj para trocar)")

    print(f"  {surface.vertex_count()} vertices, {surface.triangle_count()} triangulos, "
          f"volume {surface.volume():.6g}, estanque={surface.is_watertight()}")
    if not surface.is_watertight():
        # A surface with holes has no inside, so "carve the inside out"
        # is not a well-posed request -- better to say so than to produce
        # a mesh that is quietly wrong.
        print("  ERRO: a superficie tem furos; o mesher precisa de uma superficie fechada")
        return 1

    # 1. Volume mesh -----------------------------------------------------
    # The object becomes a hole inside a box of fluid. `margin` is how far
    # the box extends past the object, in multiples of the object's own
    # size: too small and the walls interfere with the flow, too large and
    # cells are spent on empty space.
    print("gerando malha de volume...")
    domain = aether.mesh_flow_around_object(surface, margin=2.5, background_coarsening=3.0)
    print(f"  {domain.mesh.cell_count()} celulas, "
          f"{domain.unrecovered_facet_count} facetas do objeto perdidas, "
          f"volume esculpido {domain.carved_volume:.6g} "
          f"(esperado {domain.expected_volume:.6g})")

    # 2. Boundary conditions ---------------------------------------------
    # A closed domain, driven by one box face sliding in its own plane.
    # `driving_wall_velocity` refuses a direction that is not tangential,
    # because a wall pushing along its own normal injects mass into a
    # sealed box -- unsolvable, not merely hard.
    lid = aether.driving_wall_velocity(domain, face="z_max", direction=(1.0, 0.0, 0.0), speed=1.0)

    conservation = aether.check_closed_domain_conservation(domain.mesh, lid)
    print(f"conservacao: fluxo liquido {conservation.net_flux:+.3e}, "
          f"pior face {conservation.max_face_flux:.3e} -> "
          f"{'ok' if conservation.is_conservative else 'NAO CONSERVA'}")
    if not conservation.is_conservative:
        return 1

    # 3. Will it run? ----------------------------------------------------
    # Seconds of probing instead of minutes of marching. There is no honest
    # a-priori mesh-quality threshold (see DIVIDA_TECNICA.md 4.3), so this
    # measures the two things that actually decide: that rest stays exactly
    # at rest, and the growth factor of one step.
    print("sondando estabilidade...")
    stability = aether.measure_mesh_stability(domain.mesh, viscosity=0.5)
    print(f"  raio espectral {stability.spectral_radius:.4f} "
          f"(< 1 = estavel), repouso {stability.rest_state_max_velocity:.1e}, "
          f"nao-ortogonalidade {stability.max_non_orthogonality:.2f}, "
          f"razao de volume {stability.volume_ratio:.1f}")
    if not stability.is_stable:
        print("  a sonda considera esta malha instavel; refine ou suavize a geometria")
        return 1

    # 4. Solve -----------------------------------------------------------
    solver = aether.UnstructuredCavitySolver3D(domain.mesh, 0.5, lid)
    print("marchando...")
    run = aether.run_to_steady_state(solver, domain.mesh, max_steps=600, check_every=50)
    print("  " + run.summary().replace("\n", "\n  "))

    # 5. Export ----------------------------------------------------------
    # Velocity, pressure and a derived speed, in a file ParaView or VisIt
    # opens directly. Slicing, streamlines and iso-surfaces belong there
    # rather than in this engine.
    output = "flow_around_object.vtk"
    aether.export_result_vtk(output, solver, domain.mesh)
    print(f"resultado escrito em {output} (abra no ParaView)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
