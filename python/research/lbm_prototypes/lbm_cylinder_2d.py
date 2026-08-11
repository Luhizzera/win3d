"""
Etapa 1 do simulador de tunel de vento: fluxo 2D ao redor de um cilindro,
resolvido com o Lattice Boltzmann Method (D2Q9, colisao BGK).

Serve para validar a fisica (deve aparecer o desprendimento de vortices de
Karman atras do cilindro) antes de migrar o solver para CuPy/CUDA.

Uso:
    python lbm_cylinder_2d.py
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ---------------------------------------------------------------------------
# Parametros da simulacao (tudo em unidades de lattice)
# ---------------------------------------------------------------------------
NX, NY = 400, 100                    # tamanho da grade (colunas x, linhas y)
CX, CY, R = NX // 5, NY // 2, NY // 9  # posicao/raio do cilindro
REYNOLDS = 150.0                     # numero de Reynolds alvo
U_MAX = 0.04                         # velocidade maxima de entrada (lattice)
NU = U_MAX * (2 * R) / REYNOLDS      # viscosidade cinematica (lattice)
OMEGA = 1.0 / (3.0 * NU + 0.5)       # fator de relaxacao BGK (1/tau)
N_STEPS = 8000
STEPS_PER_FRAME = 25

# ---------------------------------------------------------------------------
# Lattice D2Q9: 9 direcoes discretas de velocidade e seus pesos
# ---------------------------------------------------------------------------
C = np.array([
    [0, 0], [1, 0], [0, 1], [-1, 0], [0, -1],
    [1, 1], [-1, 1], [-1, -1], [1, -1],
])
W = np.array([4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36])
OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])          # direcao oposta (bounce-back)

RIGHT = np.array([i for i in range(9) if C[i, 0] > 0])   # direcoes que entram pela esquerda
LEFT = np.array([i for i in range(9) if C[i, 0] < 0])    # direcoes que saem pela direita
VERTICAL = np.array([i for i in range(9) if C[i, 0] == 0])

# ---------------------------------------------------------------------------
# Geometria: cilindro + paredes solidas (topo/base) do "tunel"
# ---------------------------------------------------------------------------
X, Y = np.meshgrid(np.arange(NX), np.arange(NY), indexing="ij")
cylinder = (X - CX) ** 2 + (Y - CY) ** 2 <= R ** 2
walls = (Y == 0) | (Y == NY - 1)
SOLID = cylinder | walls

# perfil de entrada tipo Poiseuille: zero nas paredes, maximo no centro
y = np.arange(NY)
U_IN = U_MAX * (1 - ((y - (NY - 1) / 2) / ((NY - 1) / 2)) ** 2)


def equilibrium(rho, ux, uy):
    cu = 3.0 * (C[:, 0, None, None] * ux + C[:, 1, None, None] * uy)
    usqr = 1.5 * (ux ** 2 + uy ** 2)
    feq = np.empty((9, NX, NY))
    for i in range(9):
        feq[i] = rho * W[i] * (1 + cu[i] + 0.5 * cu[i] ** 2 - usqr)
    return feq


def macroscopic(f):
    rho = f.sum(axis=0)
    ux = (f * C[:, 0, None, None]).sum(axis=0) / rho
    uy = (f * C[:, 1, None, None]).sum(axis=0) / rho
    return rho, ux, uy


# ---------------------------------------------------------------------------
# Inicializacao: fluido em repouso com o perfil de entrada ja aplicado
# ---------------------------------------------------------------------------
rho0 = np.ones((NX, NY))
ux0 = np.zeros((NX, NY))
uy0 = np.zeros((NX, NY))
ux0[0, :] = U_IN
f = equilibrium(rho0, ux0, uy0)


def step(f):
    # 1) saida (outflow): extrapola a ultima coluna a partir da penultima
    f[LEFT, -1, :] = f[LEFT, -2, :]

    # 2) momentos macroscopicos a partir das populacoes atuais
    rho, ux, uy = macroscopic(f)

    # 3) entrada (Zou/He): impoe velocidade e reconstroi populacoes desconhecidas
    ux[0, :] = U_IN
    uy[0, :] = 0.0
    rho[0, :] = (
        f[VERTICAL, 0, :].sum(axis=0) + 2 * f[LEFT, 0, :].sum(axis=0)
    ) / (1 - ux[0, :])

    feq = equilibrium(rho, ux, uy)
    f[RIGHT, 0, :] = feq[RIGHT, 0, :] + (f[OPP[RIGHT], 0, :] - feq[OPP[RIGHT], 0, :])

    # 4) colisao BGK
    feq = equilibrium(rho, ux, uy)
    f_post = f - OMEGA * (f - feq)

    # 5) bounce-back nas celulas solidas (cilindro + paredes), com f pre-colisao
    for i in range(9):
        f_post[i, SOLID] = f[OPP[i], SOLID]

    # 6) streaming: cada populacao se desloca na sua direcao
    f_new = np.empty_like(f_post)
    for i in range(9):
        f_new[i] = np.roll(np.roll(f_post[i], C[i, 0], axis=0), C[i, 1], axis=1)

    return f_new


def main():
    global f

    fig, ax = plt.subplots(figsize=(9, 3))
    speed0 = np.zeros((NY, NX))
    im = ax.imshow(speed0, cmap="turbo", origin="lower", vmin=0, vmax=U_MAX * 1.6)
    ax.set_title("Velocidade do fluxo ao redor do cilindro (LBM D2Q9)")
    ax.set_xticks([])
    ax.set_yticks([])
    plt.colorbar(im, ax=ax, label="|u| (unidades de lattice)")

    def animate(frame):
        global f
        for _ in range(STEPS_PER_FRAME):
            f = step(f)
        _, ux, uy = macroscopic(f)
        speed = np.sqrt(ux ** 2 + uy ** 2)
        speed[SOLID] = np.nan
        im.set_data(speed.T)
        ax.set_title(f"Velocidade do fluxo (LBM D2Q9) - passo {frame * STEPS_PER_FRAME}")
        return [im]

    n_frames = N_STEPS // STEPS_PER_FRAME
    ani = animation.FuncAnimation(fig, animate, frames=n_frames, interval=30, blit=False)
    plt.show()


if __name__ == "__main__":
    main()
