import matplotlib
matplotlib.use("Agg")
import cupy as cp
import matplotlib.pyplot as plt

import lbm_cylinder_2d_gpu as sim

sim.benchmark(n=300)

f = sim.f
for i in range(3000):
    f = sim.step(f)
    if i % 500 == 0:
        rho, ux, uy = sim.macroscopic(f)
        assert cp.isfinite(f).all(), f"instabilidade numerica no passo {i}"
        speed = cp.sqrt(ux**2 + uy**2)
        print(f"passo {i}: max|u|={float(speed.max()):.4f}  "
              f"rho min/max={float(rho.min()):.4f}/{float(rho.max()):.4f}")

rho, ux, uy = sim.macroscopic(f)
speed = cp.sqrt(ux**2 + uy**2)
speed = cp.where(sim.SOLID, cp.nan, speed)

fig, ax = plt.subplots(figsize=(10, 3.2))
im = ax.imshow(cp.asnumpy(speed).T, cmap="turbo", origin="lower")
plt.colorbar(im, ax=ax, label="|u|")
ax.set_title("Smoke test GPU - passo 3000")
fig.savefig("smoke_test_gpu_output.png", dpi=120)
print("OK - imagem salva em smoke_test_gpu_output.png")
