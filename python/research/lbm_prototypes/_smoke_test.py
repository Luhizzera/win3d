"""Teste rapido e headless do solver LBM: roda alguns passos e checa
estabilidade numerica (sem NaN/Inf), depois salva um PNG do campo de
velocidade pra inspecao visual."""

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

import lbm_cylinder_2d as sim

f = sim.f
for i in range(3000):
    f = sim.step(f)
    if i % 500 == 0:
        rho, ux, uy = sim.macroscopic(f)
        assert np.isfinite(f).all(), f"instabilidade numerica no passo {i}"
        speed = np.sqrt(ux**2 + uy**2)
        print(f"passo {i}: max|u|={speed.max():.4f}  rho min/max={rho.min():.4f}/{rho.max():.4f}")

rho, ux, uy = sim.macroscopic(f)
speed = np.sqrt(ux**2 + uy**2)
speed[sim.SOLID] = np.nan

fig, ax = plt.subplots(figsize=(9, 3))
im = ax.imshow(speed.T, cmap="turbo", origin="lower")
plt.colorbar(im, ax=ax, label="|u|")
ax.set_title("Smoke test - passo 3000")
fig.savefig("smoke_test_output.png", dpi=120)
print("OK - imagem salva em smoke_test_output.png")
