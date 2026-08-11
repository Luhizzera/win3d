"""
Etapa 2: mesmo solver LBM D2Q9 da etapa 1, agora rodando na GPU NVIDIA via
CuPy (troca quase 1:1 de numpy -> cupy). A logica fisica e identica ao
lbm_cylinder_2d.py; o que muda e o array backend (`xp`) e a necessidade de
trazer os dados de volta pra CPU (`xp.asnumpy`) so na hora de plotar.

Uso:
    python lbm_cylinder_2d_gpu.py
"""

import time

import cupy as xp  # <- unica troca estrutural em relacao a etapa 1
import numpy as np  # usado so pra montar dados que vao pro matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ---------------------------------------------------------------------------
# Parametros da simulacao (grade maior: a GPU aguenta tranquilo)
# ---------------------------------------------------------------------------
NX, NY = 800, 200
CX, CY, R = NX // 5, NY // 2, NY // 9
REYNOLDS = 150.0
U_MAX = 0.04
NU = U_MAX * (2 * R) / REYNOLDS
OMEGA = 1.0 / (3.0 * NU + 0.5)
N_STEPS = 20000
STEPS_PER_FRAME = 40

# ---------------------------------------------------------------------------
# Lattice D2Q9
# ---------------------------------------------------------------------------
C = xp.array([
    [0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
    [1, 1], [-1, 1], [-1, -1], [1, -1],
])
W = xp.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
OPP = xp.array([0, 3, 4, 1, 2, 7, 8, 5, 6])

RIGHT = xp.array([i for i in range(9) if int(C[i, 0]) > 0])
LEFT = xp.array([i for i in range(9) if int(C[i, 0]) < 0])
VERTICAL = xp.array([i for i in range(9) if int(C[i, 0]) == 0])

# ---------------------------------------------------------------------------
# Geometria
# ---------------------------------------------------------------------------
X, Y = xp.meshgrid(xp.arange(NX), xp.arange(NY), indexing="ij")
cylinder = (X - CX) ** 2 + (Y - CY) ** 2 <= R ** 2
walls = (Y == 0) | (Y == NY - 1)
SOLID = cylinder | walls

y = xp.arange(NY)
U_IN = U_MAX * (1 - ((y - (NY - 1) / 2) / ((NY - 1) / 2)) ** 2)


def equilibrium(rho, ux, uy):
    cu = 3.0 * (C[:, 0, None, None] * ux + C[:, 1, None, None] * uy)
    usqr = 1.5 * (ux ** 2 + uy ** 2)
    feq = xp.empty((9, NX, NY))
    for i in range(9):
        feq[i] = rho * W[i] * (1 + cu[i] + 0.5 * cu[i] ** 2 - usqr)
    return feq


def macroscopic(f):
    rho = f.sum(axis=0)
    ux = (f * C[:, 0, None, None]).sum(axis=0) / rho
    uy = (f * C[:, 1, None, None]).sum(axis=0) / rho
    return rho, ux, uy


rho0 = xp.ones((NX, NY))
ux0 = xp.zeros((NX, NY))
uy0 = xp.zeros((NX, NY))
ux0[0, :] = U_IN
f = equilibrium(rho0, ux0, uy0)


def step(f):
    f[LEFT, -1, :] = f[LEFT, -2, :]

    rho, ux, uy = macroscopic(f)

    ux[0, :] = U_IN
    uy[0, :] = 0.0
    rho[0, :] = (
        f[VERTICAL, 0, :].sum(axis=0) + 2 * f[LEFT, 0, :].sum(axis=0)
    ) / (1 - ux[0, :])

    feq = equilibrium(rho, ux, uy)
    f[RIGHT, 0, :] = feq[RIGHT, 0, :] + (f[OPP[RIGHT], 0, :] - feq[OPP[RIGHT], 0, :])

    feq = equilibrium(rho, ux, uy)
    f_post = f - OMEGA * (f - feq)

    for i in range(9):
        f_post[i, SOLID] = f[OPP[i], SOLID]

    f_new = xp.empty_like(f_post)
    for i in range(9):
        f_new[i] = xp.roll(xp.roll(f_post[i], int(C[i, 0]), axis=0), int(C[i, 1]), axis=1)

    return f_new


def benchmark(n=200):
    """Mede tempo por passo na GPU pra comparar com a etapa 1 (CPU)."""
    global f
    xp.cuda.Stream.null.synchronize()
    t0 = time.perf_counter()
    for _ in range(n):
        f = step(f)
    xp.cuda.Stream.null.synchronize()
    dt = time.perf_counter() - t0
    print(f"[GPU] {n} passos em {dt:.3f}s -> {dt / n * 1000:.3f} ms/passo "
          f"({NX}x{NY} celulas)")


def main():
    global f

    fig, ax = plt.subplots(figsize=(10, 3.2))
    im = ax.imshow(np.zeros((NY, NX)), cmap="turbo", origin="lower", vmin=0, vmax=U_MAX * 1.6)
    ax.set_title("Velocidade do fluxo ao redor do cilindro (LBM D2Q9, GPU/CuPy)")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, label="|u| (unidades de lattice)")

    def animate(frame):
        global f
        for _ in range(STEPS_PER_FRAME):
            f = step(f)
        _, ux, uy = macroscopic(f)
        speed = xp.sqrt(ux ** 2 + uy ** 2)
        speed = xp.where(SOLID, xp.nan, speed)
        im.set_data(xp.asnumpy(speed).T)
        ax.set_title(f"Velocidade do fluxo (LBM D2Q9, GPU) - passo {frame * STEPS_PER_FRAME}")
        return [im]

    n_frames = N_STEPS // STEPS_PER_FRAME
    ani = animation.FuncAnimation(fig, animate, frames=n_frames, interval=30, blit=False)
    plt.show()


if __name__ == "__main__":
    benchmark()
    main()
