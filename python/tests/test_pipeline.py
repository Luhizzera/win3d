"""Tests for the geometry -> volume mesh -> solver pipeline
(python/aether/pipeline.py).

Registered with ctest alongside the binding tests, for the reason
DIVIDA_TECNICA.md 5.2 gives for those: a layer nobody runs automatically is
a layer that quietly stops working. This one has a second reason -- the
pipeline is the only part of the engine that touches *imported* geometry
rather than the synthetic jittered lattice every other test builds, so it
is the only place a regression in that path would show up at all.

Run directly, it takes the extension directory and the package directory as
argv[1] and argv[2], the same convention test_unstructured_bindings.py uses.
"""

import math
import sys

failures = []


def check(condition, description):
    status = "ok " if condition else "FALHOU"
    print(f"  [{status}] {description}")
    if not condition:
        failures.append(description)


def build_icosahedron(triangle_mesh_cls, vector3_cls):
    """A unit-circumradius icosahedron: 12 vertices, 20 triangles, convex and
    watertight by construction.

    Built here rather than loaded from a checked-in STL so the test has no
    binary fixture to keep in sync, and chosen over a cube because a cube's
    square faces are coplanar quads -- exactly the facet-recovery case
    DelaunayTetrahedralization3D documents as *not* always recoverable, which
    would make a failure here ambiguous between "the pipeline broke" and
    "this geometry was never supported".
    """
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    raw = [(-1, phi, 0), (1, phi, 0), (-1, -phi, 0), (1, -phi, 0),
           (0, -1, phi), (0, 1, phi), (0, -1, -phi), (0, 1, -phi),
           (phi, 0, -1), (phi, 0, 1), (-phi, 0, -1), (-phi, 0, 1)]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11),
             (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9),
             (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    norm = math.sqrt(1.0 + phi * phi)
    mesh = triangle_mesh_cls()
    for x, y, z in raw:
        mesh.add_vertex(vector3_cls(x / norm, y / norm, z / norm))
    for a, b, c in faces:
        mesh.add_triangle(a, b, c)
    return mesh


def test_mesh_generation_is_geometrically_exact(aether, solid):
    """The fluid domain is the box minus the object, and that is an identity
    rather than an approximation: every tetrahedron's volume is exact (the
    determinant formula, no quadrature), so the carved total must match
    box-volume minus object-volume to roundoff, not to a tolerance chosen
    for convenience.
    """
    print("malha a partir de geometria importada: exatidao geometrica")
    # background_coarsening=3.0 rather than the default: the claims below
    # are all exact identities or sign checks, none of which need
    # resolution, and tetrahedralize() is still O(N^2) -- so a finer mesh
    # would multiply this suite's runtime without strengthening a single
    # assertion.
    domain = aether.mesh_flow_around_object(solid, margin=2.5, background_coarsening=3.0)

    check(domain.mesh.cell_count() > 0, f"a malha tem celulas ({domain.mesh.cell_count()})")
    check(domain.unrecovered_facet_count == 0,
          f"toda faceta do objeto foi recuperada ({domain.unrecovered_facet_count} perdidas)")
    relative = abs(domain.carved_volume - domain.expected_volume) / domain.expected_volume
    check(relative < 1e-12, f"volume esculpido = caixa - objeto (erro relativo {relative:.2e})")

    # Every one of the object's own faces must survive as a boundary face:
    # if remove_region had leaked past a wall, the object would not be a
    # hole in the domain and these would be interior faces instead.
    labels = {}
    for f in range(domain.mesh.face_count()):
        if not domain.mesh.is_boundary_face(f):
            continue
        label = aether.classify_boundary_face(domain.mesh.face(f).centroid, domain)
        labels[label] = labels.get(label, 0) + 1
    check(labels.get("object", 0) == solid.triangle_count(),
          f"as {solid.triangle_count()} faces do objeto continuam no contorno "
          f"({labels.get('object', 0)})")
    check(all(labels.get(side, 0) > 0
              for side in ("x_min", "x_max", "y_min", "y_max", "z_min", "z_max")),
          "as seis faces da caixa aparecem no contorno")
    return domain


def test_conservation_check_catches_non_tangential_lid(aether, domain):
    """The check exists because this exact mistake was made and cost real
    time: a wall velocity with a component along its own face's normal
    injects mass into a sealed box. The test does not just confirm the
    checker passes a good case -- a checker that always returns "fine"
    would do that -- it builds the bad case deliberately and requires it to
    be caught.
    """
    print("checagem de conservacao: pega tampa nao-tangencial")

    good = aether.driving_wall_velocity(domain, face="z_max", direction=(1.0, 0.0, 0.0))
    report = aether.check_closed_domain_conservation(domain.mesh, good)
    check(report.is_conservative,
          f"tampa tangencial conserva massa (desbalanco relativo {report.relative_imbalance:.2e})")
    # Exactly zero, not merely small: a tangential lid's velocity is
    # perpendicular to its own face's area vector on every single face, so
    # each dot product vanishes identically rather than cancelling in a sum.
    check(report.max_face_flux == 0.0,
          f"nenhuma parede deixa massa passar ({report.max_face_flux:.3e} no pior caso)")

    # driving_wall_velocity refuses to build the bad case at all, so it is
    # constructed by hand here -- the failure mode has to be reachable for
    # the checker to be worth anything.
    def normal_lid(position):
        if aether.classify_boundary_face(position, domain) != "x_min":
            return aether.Vector3(0.0, 0.0, 0.0)
        return aether.Vector3(1.0, 0.0, 0.0)  # normal to x_min: pushes mass in

    bad = aether.check_closed_domain_conservation(domain.mesh, normal_lid)
    check(not bad.is_conservative,
          f"tampa normal a face e reprovada (desbalanco relativo {bad.relative_imbalance:.2e})")
    check(bad.max_face_flux > 0.0,
          f"e a parede permeavel e apontada ({bad.max_face_flux:.3e} de fluxo numa face)")

    raised = False
    try:
        aether.driving_wall_velocity(domain, face="x_min", direction=(1.0, 0.0, 0.0))
    except ValueError:
        raised = True
    check(raised, "driving_wall_velocity recusa direcao nao-tangencial em vez de corrigi-la")


def test_probe_predicts_the_run_and_the_run_shows_a_vortex(aether, domain):
    """Two claims that share one march, because the march is what this file
    costs and running it twice bought nothing.

    **The probe predicts.** Its whole purpose is to say in seconds what a
    long run would reveal in minutes, so it is checked against the run: the
    probe must call this mesh stable, and the mesh must then actually be
    stable. The rest-state half is a bit-exact claim, not a tolerance --
    with every wall stationary there is no forcing anywhere, so any nonzero
    velocity is a structural defect rather than accumulated error.

    **The run shows the vortex.** The same claim every cavity solver in
    this project is validated against, now on a mesh built from imported
    geometry with an object carved out of it: fluid near the moving lid
    follows it, and fluid near the opposite wall must return against it. In
    a sealed box that reversal is forced by mass conservation, not merely
    plausible -- which is what makes it checkable without any reference
    solution.
    """
    print("sonda preve a marcha, e a marcha mostra o vortice")
    report = aether.measure_mesh_stability(domain.mesh, viscosity=0.5)

    check(report.rest_state_max_velocity == 0.0,
          f"repouso permanece exatamente zero ({report.rest_state_max_velocity:.3e})")
    check(report.spectral_radius < 1.0,
          f"raio espectral abaixo de 1 ({report.spectral_radius:.4f})")
    check(report.is_stable, "a sonda considera a malha estavel")
    print(f"       {report.cell_count} celulas, naoOrtog={report.max_non_orthogonality:.2f}, "
          f"estencilDeficiente={report.deficient_stencil_count}, "
          f"razaoVol={report.volume_ratio:.1f}")

    lid = aether.driving_wall_velocity(domain, face="z_max", direction=(1.0, 0.0, 0.0), speed=1.0)
    solver = aether.UnstructuredCavitySolver3D(domain.mesh, 0.5, lid)
    run = aether.run_to_steady_state(solver, domain.mesh, max_steps=300, check_every=50)
    print("       " + run.summary().replace("\n", "\n       "))

    check(not run.diverged, "a marcha real nao divergiu, como a sonda previu")
    # A sealed box exchanges no mass with anything, at any point in the run.
    check(abs(run.net_boundary_flux) < 1e-12,
          f"caixa selada nao troca massa ({run.net_boundary_flux:+.2e})")
    # Nothing in a closed cavity can outrun the lid that drives it.
    check(run.max_velocity < 1.0,
          f"velocidade maxima abaixo da tampa ({run.max_velocity:.4f} < 1.0)")

    span = domain.box_max.z - domain.box_min.z
    top, bottom = [], []
    for c in range(domain.mesh.cell_count()):
        z = domain.mesh.cell_centroid(c).z
        entry = (domain.mesh.cell_volume(c), solver.velocity(c).x)
        if z > domain.box_max.z - 0.2 * span:
            top.append(entry)
        elif z < domain.box_min.z + 0.2 * span:
            bottom.append(entry)

    top_mean = sum(v * u for v, u in top) / sum(v for v, _ in top)
    bottom_mean = sum(v * u for v, u in bottom) / sum(v for v, _ in bottom)
    print(f"       u medio topo={top_mean:+.5f} fundo={bottom_mean:+.5f}")
    check(top_mean > 0.0, f"o fluido junto a tampa a segue ({top_mean:+.5f})")
    check(bottom_mean < 0.0, f"o retorno junto ao fundo e contrario ({bottom_mean:+.5f})")
    return solver


def test_vtk_export_carries_the_result(aether, domain, solver):
    """The VTK round-trip itself is checked in C++ (postprocessing_tests);
    what is checked here is the layer above it -- that
    `export_result_vtk` hands the *solver's own* fields over, with the
    cell counts and the field names a post-processor will look for, rather
    than writing a well-formed file about the wrong data.
    """
    print("exportacao VTK: leva o resultado do solver")
    import os
    import tempfile

    path = os.path.join(tempfile.gettempdir(), "aether_pipeline_result.vtk")
    aether.export_result_vtk(path, solver, domain.mesh)
    check(os.path.exists(path), "o arquivo VTK foi escrito")

    with open(path, "r", encoding="ascii") as handle:
        text = handle.read()
    os.remove(path)

    check(text.startswith("# vtk DataFile Version"), "cabecalho VTK reconhecivel")
    check(f"CELLS {domain.mesh.cell_count()}" in text,
          f"declara as {domain.mesh.cell_count()} celulas da malha")
    for field in ("pressure", "speed", "velocity"):
        check(field in text, f"o campo '{field}' esta no arquivo")

    # The written values must be the ones the solver actually holds -- a
    # file with the right shape and someone else's numbers would pass every
    # check above. Comparing the whole pressure block rather than its first
    # entry, because cell 0 is the *pinned* reference cell in a closed
    # domain and therefore holds exactly 0: an assertion that only looked
    # at it would be comparing 0 against 0 and could not fail.
    marker = "LOOKUP_TABLE default\n"
    start = text.index(marker) + len(marker)
    lines = text[start:].split("\n")
    written = [float(lines[c]) for c in range(domain.mesh.cell_count())]
    expected = [solver.pressure(c) for c in range(domain.mesh.cell_count())]
    check(written == expected,
          f"as {len(written)} pressoes escritas sao exatamente as do solver")
    check(any(value != 0.0 for value in written),
          f"e o campo escrito nao e trivialmente nulo (max |p| = {max(abs(v) for v in written):.3e})")


def main():
    if len(sys.argv) > 1:
        sys.path.insert(0, sys.argv[1])
    if len(sys.argv) > 2:
        sys.path.insert(0, sys.argv[2])

    import aether

    solid = build_icosahedron(aether.TriangleMesh, aether.Vector3)
    check(solid.is_watertight(), f"o objeto de teste e estanque ({solid.triangle_count()} triangulos)")

    domain = test_mesh_generation_is_geometrically_exact(aether, solid)
    test_conservation_check_catches_non_tangential_lid(aether, domain)
    solver = test_probe_predicts_the_run_and_the_run_shows_a_vortex(aether, domain)
    test_vtk_export_carries_the_result(aether, domain, solver)

    if failures:
        print(f"\nFALHOU: {len(failures)} verificacao(oes)")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("\ntest_pipeline: todas as verificacoes passaram")
    return 0


if __name__ == "__main__":
    sys.exit(main())
