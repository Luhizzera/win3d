"""
Etapa 3: mesmo modelo fisico (LBM D2Q9, BGK) das etapas 1 e 2, mas agora com
um kernel CUDA escrito a mao (via Numba) que funde num unico lancamento por
passo: colisao + condicoes de contorno (entrada/saida/paredes) + streaming.

Isso elimina o overhead das ~30 chamadas de kernel pequenas que a versao
CuPy fazia por passo (uma por operacao de array), que era o gargalo
identificado na etapa 2.

Estrategia do kernel (push/scatter com buffer duplo):
    - cada thread cuida de UMA celula da grade
    - le suas 9 populacoes de f_in
    - se for celula solida: reflete (bounce-back) e "empurra" o resultado
      para as celulas vizinhas em f_out
    - se for celula de entrada/saida: corrige as populacoes desconhecidas
      (Zou/He na entrada, extrapolacao na saida)
    - faz a colisao BGK e empurra cada uma das 9 populacoes resultantes
      para a celula vizinha correspondente em f_out
    Como cada celula de destino em f_out[i] so pode receber de UMA celula de
    origem por direcao, nao ha condicao de corrida entre threads.

Uso:
    python lbm_cylinder_2d_cuda.py
"""

import os
import sys
import time

# aponta o numba pros pacotes de CUDA toolkit instalados via pip (nvvm +
# cudart), ja que nao estamos usando um ambiente conda nem um toolkit
# instalado no sistema. Precisa ser feito ANTES de importar numba.cuda.
_VENV_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(sys.executable)))
_CUDA_HOME = os.path.join(_VENV_ROOT, "cuda_home")
if os.path.isdir(_CUDA_HOME):
    os.environ.setdefault("CUDA_HOME", _CUDA_HOME)

import numpy as np
import numba
from numba import cuda
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ---------------------------------------------------------------------------
# Parametros da simulacao
# ---------------------------------------------------------------------------
NX, NY = 800, 200
CX_OBJ, CY_OBJ, R_OBJ = NX // 5, NY // 2, NY // 9
REYNOLDS = 150.0
U_MAX = 0.04
NU = U_MAX * (2 * R_OBJ) / REYNOLDS
OMEGA = 1.0 / (3.0 * NU + 0.5)
N_STEPS = 30000
STEPS_PER_FRAME = 60

# ---------------------------------------------------------------------------
# Lattice D2Q9 (constantes compile-time, usadas dentro do kernel)
# ---------------------------------------------------------------------------
CX = (0, 1, 0, -1, 0, 1, -1, -1, 1)
CY = (0, 0, 1, 0, -1, 1, 1, -1, -1)
WT = (4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36)
OPP = (0, 3, 4, 1, 2, 7, 8, 5, 6)
RIGHT = (1, 5, 8)   # c_x > 0: desconhecidas na entrada
LEFT = (3, 6, 7)    # c_x < 0: desconhecidas na saida
VERTICAL = (0, 2, 4)  # c_x == 0

THREADS = (16, 16)
BLOCKS = (
    (NX + THREADS[0] - 1) // THREADS[0],
    (NY + THREADS[1] - 1) // THREADS[1],
)


@cuda.jit
def lbm_step_kernel(f_in, f_out, solid, u_in, nx, ny, omega):
    x, y = cuda.grid(2)
    if x >= nx or y >= ny:
        return

    lf = cuda.local.array(9, numba.float64)
    for i in range(9):
        lf[i] = f_in[i, x, y]

    # --- celula solida: bounce-back puro, sem colisao ----------------------
    if solid[x, y]:
        for i in range(9):
            xf = x + CX[i]
            yf = y + CY[i]
            if 0 <= xf < nx and 0 <= yf < ny:
                f_out[i, xf, yf] = lf[OPP[i]]
        return

    # --- momentos macroscopicos a partir das populacoes recebidas ---------
    rho = 0.0
    ux = 0.0
    uy = 0.0
    for i in range(9):
        rho += lf[i]
        ux += lf[i] * CX[i]
        uy += lf[i] * CY[i]
    ux /= rho
    uy /= rho

    # --- contorno de entrada (Zou/He): impoe velocidade --------------------
    if x == 0:
        ux = u_in[y]
        uy = 0.0
        known = 0.0
        for i in VERTICAL:
            known += lf[i]
        for i in LEFT:
            known += 2.0 * lf[i]
        rho = known / (1.0 - ux)
        usqr = 1.5 * (ux * ux + uy * uy)
        for i in RIGHT:
            j = OPP[i]
            cu_i = 3.0 * (CX[i] * ux + CY[i] * uy)
            feq_i = rho * WT[i] * (1.0 + cu_i + 0.5 * cu_i * cu_i - usqr)
            cu_j = 3.0 * (CX[j] * ux + CY[j] * uy)
            feq_j = rho * WT[j] * (1.0 + cu_j + 0.5 * cu_j * cu_j - usqr)
            lf[i] = feq_i + (lf[j] - feq_j)

    # --- contorno de saida: extrapolacao (gradiente zero) ------------------
    elif x == nx - 1:
        for i in LEFT:
            lf[i] = f_in[i, nx - 2, y]
        rho = 0.0
        ux = 0.0
        uy = 0.0
        for i in range(9):
            rho += lf[i]
            ux += lf[i] * CX[i]
            uy += lf[i] * CY[i]
        ux /= rho
        uy /= rho

    # --- colisao BGK + streaming (push) ------------------------------------
    usqr = 1.5 * (ux * ux + uy * uy)
    for i in range(9):
        cu = 3.0 * (CX[i] * ux + CY[i] * uy)
        feq = rho * WT[i] * (1.0 + cu + 0.5 * cu * cu - usqr)
        f_star = lf[i] - omega * (lf[i] - feq)
        xf = x + CX[i]
        yf = y + CY[i]
        if 0 <= xf < nx and 0 <= yf < ny:
            f_out[i, xf, yf] = f_star


def make_initial_state():
    X, Y = np.meshgrid(np.arange(NX), np.arange(NY), indexing="ij")
    cylinder = (X - CX_OBJ) ** 2 + (Y - CY_OBJ) ** 2 <= R_OBJ ** 2
    walls = (Y == 0) | (Y == NY - 1)
    solid = cylinder | walls

    y = np.arange(NY)
    u_in = U_MAX * (1 - ((y - (NY - 1) / 2) / ((NY - 1) / 2)) ** 2)

    rho0 = np.ones((NX, NY))
    ux0 = np.zeros((NX, NY))
    uy0 = np.zeros((NX, NY))
    ux0[0, :] = u_in

    cu = 3.0 * (np.array(CX)[:, None, None] * ux0 + np.array(CY)[:, None, None] * uy0)
    usqr = 1.5 * (ux0 ** 2 + uy0 ** 2)
    f0 = np.empty((9, NX, NY))
    for i in range(9):
        f0[i] = rho0 * WT[i] * (1 + cu[i] + 0.5 * cu[i] ** 2 - usqr)

    return f0, solid, u_in


class Simulation:
    def __init__(self):
        f0, solid, u_in = make_initial_state()
        self.d_f_in = cuda.to_device(np.ascontiguousarray(f0))
        self.d_f_out = cuda.to_device(np.zeros_like(f0))
        self.d_solid = cuda.to_device(np.ascontiguousarray(solid))
        self.d_u_in = cuda.to_device(np.ascontiguousarray(u_in))
        self.solid_host = solid

    def step(self):
        lbm_step_kernel[BLOCKS, THREADS](
            self.d_f_in, self.d_f_out, self.d_solid, self.d_u_in, NX, NY, OMEGA
        )
        self.d_f_in, self.d_f_out = self.d_f_out, self.d_f_in

    def speed_field(self):
        f = self.d_f_in.copy_to_host()
        rho = f.sum(axis=0)
        ux = (f * np.array(CX)[:, None, None]).sum(axis=0) / rho
        uy = (f * np.array(CY)[:, None, None]).sum(axis=0) / rho

        # As colunas x=0 (entrada) e x=NX-1 (saida) tem 3 dos 9 canais
        # deliberadamente vazios no buffer bruto: o kernel reconstroi essas
        # populacoes localmente (Zou/He / extrapolacao) e as empurra direto
        # para a coluna vizinha, sem jamais armazena-las de volta na coluna
        # de contorno. Isso e correto para a fisica (a coluna vizinha recebe
        # o valor certo), mas significa que somar o buffer bruto nessas duas
        # colunas especificas subestima rho/u. Para fins de leitura/plot,
        # copiamos da coluna interna adjacente.
        rho[0, :] = rho[1, :]
        ux[0, :] = ux[1, :]
        uy[0, :] = uy[1, :]
        rho[-1, :] = rho[-2, :]
        ux[-1, :] = ux[-2, :]
        uy[-1, :] = uy[-2, :]

        speed = np.sqrt(ux ** 2 + uy ** 2)
        return rho, speed


def benchmark(sim, n=300):
    cuda.synchronize()
    t0 = time.perf_counter()
    for _ in range(n):
        sim.step()
    cuda.synchronize()
    dt = time.perf_counter() - t0
    print(f"[CUDA kernel] {n} passos em {dt:.3f}s -> {dt / n * 1000:.3f} ms/passo "
          f"({NX}x{NY} celulas)")


def main():
    sim = Simulation()

    fig, ax = plt.subplots(figsize=(10, 3.2))
    im = ax.imshow(np.zeros((NY, NX)), cmap="turbo", origin="lower", vmin=0, vmax=U_MAX * 1.6)
    ax.set_title("Velocidade do fluxo (LBM D2Q9, kernel CUDA customizado)")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, label="|u| (unidades de lattice)")

    def animate(frame):
        for _ in range(STEPS_PER_FRAME):
            sim.step()
        _, speed = sim.speed_field()
        speed[sim.solid_host] = np.nan
        im.set_data(speed.T)
        ax.set_title(f"Velocidade do fluxo (kernel CUDA) - passo {frame * STEPS_PER_FRAME}")
        return [im]

    n_frames = N_STEPS // STEPS_PER_FRAME
    ani = animation.FuncAnimation(fig, animate, frames=n_frames, interval=30, blit=False)
    plt.show()


if __name__ == "__main__":
    sim = Simulation()
    benchmark(sim, n=300)
    main()
