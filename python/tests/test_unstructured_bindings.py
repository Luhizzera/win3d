"""Python-level gate for the unstructured FVM layer (DIVIDA_TECNICA.md 5.2).

**Why this file exists, and why it is a ctest and not a notebook.** Bindings
that nothing runs are the same failure mode as a CI workflow that never
executed: a hypothesis, not a safety net. Every other layer of this engine is
exercised from C++ by ctest, and the binding modules were until now only ever
*compiled* -- which catches a signature error and nothing else. A binding can
compile and still hand back a dangling pointer, drop a default argument, or
mis-order constructor parameters, and all three are mistakes this file would
catch.

**What it checks, and what it deliberately does not.** It does not re-verify
the physics: the Fourier-series comparison and the mesh convergence study
belong in engine/solver/tests, where they already are, and duplicating them
here would mean two things to keep in step. What it checks is the property
the C++ suite structurally cannot: that the *binding surface* works --
callables crossing into C++, ownership surviving Python's garbage collector,
diagnostics reaching Python with the values C++ computed.

The physical claims it does make are the ones that are free: consequences of
the equations that any correct solver must satisfy and that need no reference
table recalled from memory. The maximum principle for the Laplacian, and the
mass balance of a channel with an outlet.

Run: python test_unstructured_bindings.py <extension-dir> <python-package-dir>
"""

import gc
import math
import random
import sys

failures = []


def check(condition, description):
    """Report and record rather than raising, so one run reports every
    failure instead of only the first -- a binding break tends to come in
    groups, and finding them one build at a time is the slow way."""
    status = "ok " if condition else "FALHOU"
    print(f"  [{status}] {description}")
    if not condition:
        failures.append(description)


def build_jittered_lattice(n, seed=17):
    """A cube lattice with the interior points perturbed, tetrahedralized.

    Jittered rather than regular on purpose: a regular lattice tetrahedralizes
    into a mesh far more orthogonal than anything a real geometry produces, so
    it would flatter every gradient scheme here. Boundary points are left
    exactly on the faces so the domain stays the unit cube -- the boundary
    selectors below identify faces by plane, and a jittered boundary point
    would put a face centroid off that plane.
    """
    import aether_mesh_py as mesh_module

    tets = mesh_module.DelaunayTetrahedralization3D()
    rng = random.Random(seed)
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(n + 1):
                point = []
                for index in (i, j, k):
                    value = index / n
                    if 0 < index < n:
                        value += rng.uniform(-0.25, 0.25) / n
                    point.append(value)
                tets.add_point(*point)
    tets.tetrahedralize()
    return mesh_module.TetrahedralMesh.from_tetrahedralization(tets)


def test_mesh_surface(mesh):
    print("mesh: identidades geometricas")
    check(mesh.cell_count() > 0, f"malha tem celulas ({mesh.cell_count()})")
    check(abs(mesh.total_volume() - 1.0) < 1e-9,
          f"volume total = 1 no cubo unitario ({mesh.total_volume():.12f})")
    check(mesh.boundary_face_count() > 0,
          f"tem faces de contorno ({mesh.boundary_face_count()})")

    # The discrete divergence theorem: a closed polyhedron's outward area
    # vectors sum to zero identically, not approximately. A nonzero result
    # means a normal is misoriented or a face is missing -- the two bugs the
    # connectivity could plausibly have.
    worst = 0.0
    for cell in range(mesh.cell_count()):
        total = mesh.cell_area_vector_sum(cell)
        worst = max(worst, math.sqrt(total.x**2 + total.y**2 + total.z**2))
    check(worst < 1e-12, f"soma dos vetores de area = 0 por celula (pior {worst:.2e})")

    # Face handles cross the boundary intact, including the sentinel.
    boundary = next(f for f in range(mesh.face_count()) if mesh.is_boundary_face(f))
    interior = next(f for f in range(mesh.face_count()) if not mesh.is_boundary_face(f))
    from aether_mesh_py import TetrahedralMesh

    check(mesh.face(boundary).neighbour == TetrahedralMesh.NO_NEIGHBOUR,
          "face de contorno traz o sentinela NO_NEIGHBOUR")
    check(mesh.face(interior).neighbour != TetrahedralMesh.NO_NEIGHBOUR,
          "face interior traz um vizinho de verdade")
    check(len(mesh.face(interior).vertices) == 3, "face tem tres vertices")


def test_diffusion_maximum_principle(mesh):
    """Laplace's equation obeys a maximum principle: with no source term the
    solution is bounded by its own boundary data, everywhere. Nothing about
    the mesh or the scheme changes that, so it holds at any resolution and
    needs no reference solution -- which is exactly the kind of claim this
    project prefers to check.
    """
    import aether_solver_py as solver_module

    print("difusao: principio do maximo e superficie de diagnostico")
    hot, cold = 100.0, 0.0
    solver = solver_module.UnstructuredDiffusionSolver(mesh)
    solver.set_dirichlet_boundary(lambda c: c.y > 1.0 - 1e-9, hot)
    solver.set_dirichlet_boundary(lambda c: c.y < 1e-9, cold)
    solver.set_dirichlet_boundary(lambda c: c.x < 1e-9, cold)
    solver.set_dirichlet_boundary(lambda c: c.x > 1.0 - 1e-9, cold)

    sweeps = solver.solve_conjugate_gradient()
    values = [solver.value(cell) for cell in range(solver.cell_count())]

    check(solver.cell_count() == mesh.cell_count(), "cell_count concorda com a malha")
    check(sweeps > 0, f"rodou pelo menos uma varredura externa ({sweeps})")
    check(min(values) >= cold - 1e-6,
          f"nenhuma celula abaixo do contorno frio (min {min(values):.6f})")
    check(max(values) <= hot + 1e-6,
          f"nenhuma celula acima do contorno quente (max {max(values):.6f})")
    check(max(values) > cold + 1.0, f"o campo nao e trivial (max {max(values):.3f})")

    # The hot edge is at y = 1, so the field has to increase with y on
    # average. Comparing halves rather than pointwise: individual tetrahedra
    # are irregular, the claim is about the field.
    upper = [v for cell, v in enumerate(values) if mesh.cell_centroid(cell).y > 0.5]
    lower = [v for cell, v in enumerate(values) if mesh.cell_centroid(cell).y <= 0.5]
    check(sum(upper) / len(upper) > sum(lower) / len(lower),
          f"metade quente mais quente ({sum(upper)/len(upper):.3f} > {sum(lower)/len(lower):.3f})")

    gradients = solver.cell_gradients()
    check(len(gradients) == mesh.cell_count(), "cell_gradients traz um vetor por celula")
    check(any(abs(g.y) > 1e-6 for g in gradients), "o gradiente em y nao e nulo")

    # Diagnostics have to arrive with real values, not defaults: a binding
    # that returned 0.0 for every one of these would look like a clean mesh.
    check(solver.max_non_orthogonality() > 0.0,
          f"nao-ortogonalidade medida ({solver.max_non_orthogonality():.3f})")
    print(f"       {solver.cell_count()} celulas, {sweeps} varreduras, "
          f"naoOrtog={solver.max_non_orthogonality():.3f}, "
          f"estencilDeficiente={solver.deficient_stencil_count()}, "
          f"mudancaExterna={solver.last_outer_change():.2e}")


def test_cavity_topology(mesh):
    """The same claim the C++ suite makes for this solver, restated through
    the bindings: at low Reynolds number a lid-driven cavity has one primary
    vortex, so fluid near the lid follows it and fluid near the floor must
    move against it. Mass cannot pile up in a closed box, so this is forced
    rather than merely plausible.
    """
    import aether_core_py as core_module
    import aether_solver_py as solver_module

    print("cavidade fechada: topologia do vortice primario")
    lid_speed = 1.0

    def wall_velocity(p):
        if p.z < 1.0 - 1e-9:
            return core_module.Vector3(0.0, 0.0, 0.0)
        # Tapered to zero at the edges it meets, which removes the corner
        # discontinuity a constant lid would impose.
        taper = math.sin(math.pi * p.x) ** 2 * math.sin(math.pi * p.y) ** 2
        return core_module.Vector3(lid_speed * taper, 0.0, 0.0)

    # is_outlet omitted entirely: exercises the default argument, and a
    # closed cavity is the case where the solver has to pin a reference cell
    # instead of leaning on a prescribed outlet pressure.
    solver = solver_module.UnstructuredCavitySolver3D(mesh, 0.1, wall_velocity)

    dt = solver.stable_time_step()
    check(dt > 0.0, f"passo estavel positivo ({dt:.2e})")
    steps = int(2.0 / dt)
    for _ in range(steps):
        solver.step(dt)

    check(abs(solver.time() - steps * dt) < 1e-9 * max(1.0, steps * dt),
          "o tempo simulado acompanha os passos dados")

    top = [(mesh.cell_volume(c), solver.velocity(c).x) for c in range(mesh.cell_count())
           if mesh.cell_centroid(c).z > 0.8]
    bottom = [(mesh.cell_volume(c), solver.velocity(c).x) for c in range(mesh.cell_count())
              if mesh.cell_centroid(c).z < 0.2]
    top_mean = sum(v * u for v, u in top) / sum(v for v, _ in top)
    bottom_mean = sum(v * u for v, u in bottom) / sum(v for v, _ in bottom)

    print(f"       {steps} passos ate t={solver.time():.3f}, u medio topo={top_mean:+.5f} "
          f"fundo={bottom_mean:+.5f}, divergencia={solver.max_face_divergence():.3e}, "
          f"estencilDeficiente={solver.deficient_stencil_count()}")
    check(top_mean > 0.01, f"o fluido junto a tampa a segue ({top_mean:+.5f})")
    check(bottom_mean < 0.0, f"o retorno junto ao fundo e contrario ({bottom_mean:+.5f})")
    # A sealed box can gain or lose nothing at all.
    check(abs(solver.net_boundary_flux()) < 1e-12,
          f"caixa fechada nao troca massa ({solver.net_boundary_flux():+.2e})")


def run_channel(mesh, correctors, end_time=2.0):
    import aether_core_py as core_module
    import aether_solver_py as solver_module

    inlet_speed = 1.0

    def wall_velocity(p):
        if p.x < 1e-9:
            return core_module.Vector3(inlet_speed, 0.0, 0.0)
        return core_module.Vector3(0.0, 0.0, 0.0)

    solver = solver_module.UnstructuredCavitySolver3D(
        mesh, 0.1, wall_velocity, lambda p: p.x > 1.0 - 1e-9, 0.0,
        pressure_correctors=correctors)

    dt = solver.stable_time_step()
    for _ in range(int(end_time / dt)):
        solver.step(dt)

    inflow = sum(abs(mesh.face(f).area_vector.x) * inlet_speed
                 for f in range(mesh.face_count())
                 if mesh.is_boundary_face(f) and mesh.face(f).centroid.x < 1e-9)
    return solver, inflow


def test_channel_mass_balance(mesh):
    """Whatever enters must leave. A global statement no scheme satisfies by
    accident, and the one a sealed domain cannot satisfy at all -- so it is
    also the check that the outlet is not silently still a wall.

    **And it is the number that pays for these bindings.** The first run of
    this file failed here: the default four non-orthogonal correctors left
    2.7e-04 of the inflow unaccounted on a jittered mesh, where the C++
    suite's own channel case closes to 1e-13. Diagnosing that took three
    Python scripts and no rebuilds -- ruling out dropped faces (there are
    none), ruling out an undecayed transient (t = 2, 4 and 8 all sit at
    1e-04), and finally sweeping the corrector count, which is what showed
    the residual halving per corrector. That investigation is exactly the one
    DIVIDA_TECNICA.md 5.2 said the bindings existed to make cheap.

    So this checks the *relationship* rather than a single number: more
    correctors must mean a smaller imbalance, and enough of them must reach
    machine precision. A fixed threshold at the default count would have
    encoded the mesh this test happens to build.
    """
    print("canal com saida: balanco global de massa vs. correctores")

    results = []
    for correctors in (4, 16, 64):
        solver, inflow = run_channel(mesh, correctors)
        net = solver.net_boundary_flux()
        results.append((correctors, abs(net), solver))
        print(f"       {correctors:2d} correctores: balanco={net:+.3e} "
              f"({100.0 * abs(net) / inflow:.4f}% da entrada), "
              f"mudanca de pressao={solver.last_pressure_change():.2e}")

    check(inflow > 0.0, f"a entrada esta de fato injetando massa ({inflow:.5f})")

    # Monotone in the corrector count: this is the claim that the deferred
    # correction is converging rather than fighting something, and it is what
    # makes the count a knob instead of a mystery.
    check(results[1][1] < results[0][1],
          f"16 correctores fecham melhor que 4 ({results[1][1]:.2e} < {results[0][1]:.2e})")
    check(results[2][1] < results[1][1],
          f"64 correctores fecham melhor que 16 ({results[2][1]:.2e} < {results[1][1]:.2e})")
    check(results[2][1] < 1e-9 * inflow,
          f"com correctores suficientes fecha em precisao de maquina ({results[2][1]:.2e})")

    # The diagnostic has to agree with the outcome, or it is not a diagnostic:
    # the run that closed the balance must also report the smaller residual
    # pressure change.
    check(results[2][2].last_pressure_change() < results[0][2].last_pressure_change(),
          "last_pressure_change acompanha o balanco que ele deve prever")


def test_mesh_outlives_python_reference():
    """The solvers keep a raw pointer to the mesh, so the binding declares
    keep_alive. Without it this function is a use-after-free, and dropping
    the mesh after building the solver is the natural thing for a script to
    do -- which makes this the single most likely way these bindings could be
    wrong in a way that compiles and usually appears to work.
    """
    import aether_core_py as core_module
    import aether_solver_py as solver_module

    print("propriedade: a malha sobrevive ao seu ultimo nome em Python")
    mesh = build_jittered_lattice(2, seed=5)
    expected_cells = mesh.cell_count()
    solver = solver_module.UnstructuredCavitySolver3D(
        mesh, 0.1, lambda p: core_module.Vector3(0.0, 0.0, 0.0))

    del mesh
    gc.collect()

    # Every one of these reaches through the pointer the solver kept.
    solver.step(solver.stable_time_step())
    check(solver.cell_count() == expected_cells,
          f"o solver segue funcionando sem referencia Python a malha ({expected_cells} celulas)")
    check(math.isfinite(solver.max_face_divergence()),
          "os diagnosticos seguem finitos depois da coleta")


def main():
    if len(sys.argv) < 3:
        print("uso: test_unstructured_bindings.py <dir-das-extensoes> <dir-do-pacote>")
        return 2
    sys.path.insert(0, sys.argv[1])
    sys.path.insert(0, sys.argv[2])

    # The package's own import list is part of the surface under test: a
    # binding that exists but was never exported is not reachable by the
    # people this layer is for.
    import aether

    for name in ("TetrahedralMesh", "UnstructuredDiffusionSolver", "UnstructuredCavitySolver3D"):
        check(hasattr(aether, name), f"aether.{name} exportado pelo pacote")

    mesh = build_jittered_lattice(3)
    test_mesh_surface(mesh)
    test_diffusion_maximum_principle(mesh)
    test_cavity_topology(mesh)
    test_channel_mass_balance(mesh)
    test_mesh_outlives_python_reference()

    print()
    if failures:
        print(f"FALHOU: {len(failures)} verificacao(oes)")
        for description in failures:
            print(f"  - {description}")
        return 1
    print("todas as verificacoes passaram")
    return 0


if __name__ == "__main__":
    sys.exit(main())
