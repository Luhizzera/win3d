import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt

import lbm_cylinder_2d_cuda as sim_mod

sim = sim_mod.Simulation()
sim_mod.benchmark(sim, n=300)

for i in range(3000):
    sim.step()
    if i % 500 == 0:
        rho, speed = sim.speed_field()
        assert np.isfinite(rho).all() and np.isfinite(speed).all(), f"instabilidade no passo {i}"
        fluid_speed = speed[~sim.solid_host]
        fluid_rho = rho[~sim.solid_host]
        print(f"passo {i}: max|u|={fluid_speed.max():.4f}  rho min/max={fluid_rho.min():.4f}/{fluid_rho.max():.4f}")

rho, speed = sim.speed_field()
speed[sim.solid_host] = np.nan

fig, ax = plt.subplots(figsize=(10, 3.2))
im = ax.imshow(speed.T, cmap="turbo", origin="lower")
plt.colorbar(im, ax=ax, label="|u|")
ax.set_title("Smoke test CUDA kernel - passo 3000")
fig.savefig("smoke_test_cuda_output.png", dpi=120)
print("OK - imagem salva em smoke_test_cuda_output.png")
