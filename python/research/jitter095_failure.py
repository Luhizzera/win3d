"""What still fails at jitter 0.95, given that rest-state rho = 0.58.

CONCLUSION (2026-08-25): the remaining failure is NONLINEAR, and therefore
a different phenomenon from the linear instability item 4.3 was about and
closed. Measured: not the time step (dt x0.05 does not save it), not the
viscosity (nu x30 does not save it), but sharply amplitude-dependent --
inlet 0.01 survives 400 steps while 0.05 dies at 50. A linear instability
grows at the same rate at any amplitude; this one has a threshold, which
is the signature of a subcritical nonlinear instability.

DIVIDA_TECNICA.md 4.3 closed the linear instability around rest for outlet
domains, and the jitter sweep confirmed it: every mesh up to 0.95 has a
rest-state spectral radius below 1. Yet a *driven* run at 0.95 dies at
step 41. This decomposes that, using the instruments the item itself
built: stepWith() to switch pieces off, and loadState() to linearize
somewhere other than rest.
"""
import math
import random
import sys

import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
import aether


def lattice(n=3, seed=17, jitter=0.95):
    tets = aether.DelaunayTetrahedralization3D()
    rng = random.Random(seed)
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(n + 1):
                p = []
                for index in (i, j, k):
                    value = index / n
                    if 0 < index < n:
                        value += rng.uniform(-jitter, jitter) / n
                    p.append(value)
                tets.add_point(*p)
    tets.tetrahedralize()
    return aether.TetrahedralMesh.from_tetrahedralization(tets)


def channel(mesh, inlet=1.0, viscosity=0.1):
    return aether.UnstructuredCavitySolver3D(
        mesh, viscosity,
        lambda p: aether.Vector3(inlet, 0.0, 0.0) if p.x < 1e-9 else aether.Vector3(0.0, 0.0, 0.0),
        lambda p: p.x > 1.0 - 1e-9, 0.0)


mesh = lattice()
n = mesh.cell_count()
probe = aether.UnstructuredDiffusionSolver(mesh)
probe.set_dirichlet_boundary(lambda p: True, 0.0)
print(f"malha jitter 0.95: {n} celulas, naoOrtog={probe.max_non_orthogonality():.2f}, "
      f"estencilDeficiente={probe.deficient_stencil_count()}", flush=True)

# 1) Does a smaller dt help? The guard's own message says no; verified here
#    rather than trusted, since that claim predates the 4.3 fix.
print("\\n1) sensibilidade ao passo de tempo", flush=True)
for factor in (1.0, 0.25, 0.05):
    s = channel(mesh)
    dt = s.stable_time_step() * factor
    survived = 0
    try:
        for _ in range(400):
            s.step(dt)
            survived += 1
    except RuntimeError:
        pass
    print(f"   dt x{factor:<6}: sobreviveu {survived:4d}/400 passos"
          f"{'  (t=' + format(s.time(), '.4f') + ')' if survived else ''}", flush=True)

# 2) Which piece of the step is responsible? stepWith() exists for exactly
#    this, and it is what isolated the projection back in the 4.3 work.
print("\\n2) decomposicao do passo (stepWith)", flush=True)
for label, kw in (("completo", {}),
                  ("sem conveccao", {"convection": False}),
                  ("sem viscosidade", {"viscous": False}),
                  ("sem projecao", {"projection": False})):
    s = channel(mesh)
    dt = s.stable_time_step()
    survived = 0
    try:
        for _ in range(400):
            s.step_with(dt, **kw)
            survived += 1
    except RuntimeError:
        pass
    print(f"   {label:<18}: sobreviveu {survived:4d}/400", flush=True)

# 3) Does the driving amplitude matter? A *linear* instability grows at the
#    same rate regardless of amplitude; a nonlinear one does not. This is
#    the same amplitude test that first established the 4.3 instability was
#    linear, run again now that the linear part is fixed.
print("\\n3) sensibilidade a amplitude (linear x nao-linear)", flush=True)
for inlet in (1.0, 1e-3, 1e-6):
    s = channel(mesh, inlet=inlet)
    dt = s.stable_time_step()
    survived = 0
    try:
        for _ in range(400):
            s.step(dt)
            survived += 1
    except RuntimeError:
        pass
    print(f"   entrada {inlet:<8g}: sobreviveu {survived:4d}/400", flush=True)

# 4) Growth around the *developed* state rather than rest. If the operator
#    linearized there amplifies while the one at rest does not, the failure
#    lives in the flow the case develops, not in the mesh alone.
print("\\n4) crescimento em torno do estado desenvolvido", flush=True)
s = channel(mesh)
dt = s.stable_time_step()
steps = 0
try:
    for _ in range(30):
        s.step(dt)
        steps += 1
except RuntimeError:
    pass
if steps == 30:
    base_v = [s.velocity(c) for c in range(n)]
    base_p = [s.pressure(c) for c in range(n)]
    rng = random.Random(5)
    scale = 1e-6
    pert = [aether.Vector3(v.x + rng.uniform(-scale, scale),
                            v.y + rng.uniform(-scale, scale),
                            v.z + rng.uniform(-scale, scale)) for v in base_v]
    a = channel(mesh); a.load_state(base_v, base_p, 0.0)
    b = channel(mesh); b.load_state(pert, list(base_p), 0.0)
    prev = scale
    try:
        for i in range(12):
            a.step(dt); b.step(dt)
            diff = math.sqrt(sum((b.velocity(c).x - a.velocity(c).x) ** 2 +
                                  (b.velocity(c).y - a.velocity(c).y) ** 2 +
                                  (b.velocity(c).z - a.velocity(c).z) ** 2 for c in range(n)))
            if i >= 8:
                print(f"   passo {i:2d}: |perturbacao| = {diff:.4e}  (razao {diff/prev:.4f})", flush=True)
            prev = max(diff, 1e-300)
    except RuntimeError:
        print("   a marcha perturbada levantou antes de completar", flush=True)
else:
    print(f"   o caso base ja morre no passo {steps}, cedo demais para linearizar", flush=True)


# --- viscosity and the amplitude threshold ---------------------------

vols = sorted(mesh.cell_volume(c) for c in range(mesh.cell_count()))
print(f"razao de volume da malha: {vols[-1]/vols[0]:.1f}  (min {vols[0]:.2e}, max {vols[-1]:.2e})", flush=True)
print("\nviscosidade (Re de celula cai ao subir nu):", flush=True)
for nu in (0.1, 0.3, 1.0, 3.0):
    s = channel(mesh, viscosity=nu); dt = s.stable_time_step(); survived=0
    try:
        for _ in range(400):
            s.step(dt); survived+=1
    except RuntimeError: pass
    tag = "sobrevive" if survived==400 else f"morre em {survived}"
    print(f"   nu={nu:<5}: {tag}", flush=True)

print("\namplitude fina (onde e o limiar):", flush=True)
for inlet in (0.01, 0.05, 0.1, 0.3, 0.6):
    s = channel(mesh, inlet=inlet); dt = s.stable_time_step(); survived=0
    try:
        for _ in range(400):
            s.step(dt); survived+=1
    except RuntimeError: pass
    tag = "sobrevive" if survived==400 else f"morre em {survived}"
    print(f"   entrada {inlet:<6}: {tag}", flush=True)
