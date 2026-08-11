# LBM prototypes (experimental, not the engine's official solver)

Três estágios progressivos de um solver Lattice Boltzmann (D2Q9) para
escoamento 2D ao redor de um cilindro, usados para validar a física
(vórtices de Kármán) e explorar aceleração por GPU antes de decidir a
arquitetura definitiva:

- `lbm_cylinder_2d.py` - referência CPU/NumPy
- `lbm_cylinder_2d_gpu.py` - port para GPU via CuPy
- `lbm_cylinder_2d_cuda.py` - kernel CUDA escrito à mão (Numba), fundindo
  colisão + condições de contorno + streaming num único kernel por passo

O solver oficial da engine segue Volumes Finitos (FVM), não LBM - estes
scripts ficam aqui como referência/experimento isolado, sem dependência
com `engine/core`.
