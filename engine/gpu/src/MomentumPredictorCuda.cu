#include "aether/gpu/MomentumPredictorCuda.hpp"

#include <cuda_runtime.h>

#include <stdexcept>
#include <string>

namespace aether::gpu {

namespace {

// Staggered/cell-centered index layouts, identical to
// StaggeredCavityBase3D::indexU/indexV/indexW/indexP.
__device__ __forceinline__ std::size_t deviceIndexU(std::size_t i, std::size_t j, std::size_t k, std::size_t nx,
                                                      std::size_t ny) {
    return i + j * (nx + 1) + k * (nx + 1) * ny;
}
__device__ __forceinline__ std::size_t deviceIndexV(std::size_t i, std::size_t j, std::size_t k, std::size_t nx,
                                                      std::size_t ny) {
    return i + j * nx + k * nx * (ny + 1);
}
__device__ __forceinline__ std::size_t deviceIndexW(std::size_t i, std::size_t j, std::size_t k, std::size_t nx,
                                                      std::size_t ny) {
    return i + j * nx + k * nx * ny;
}
__device__ __forceinline__ std::size_t deviceIndexP(std::size_t i, std::size_t j, std::size_t k, std::size_t nx,
                                                      std::size_t ny) {
    return i + j * nx + k * nx * ny;
}

// Same formula as StaggeredCavityBase3D::lidVelocityAt: the lid moves in
// the u (x) direction only, tapered as sin^2 in both x and y so it
// vanishes at the walls it meets.
__device__ __forceinline__ double deviceLidVelocityAt(double x, double y, double lengthX, double lengthY,
                                                        double lidVelocity) {
    constexpr double kPi = 3.14159265358979323846;
    const double sx = sin(kPi * x / lengthX);
    const double sy = sin(kPi * y / lengthY);
    return lidVelocity * sx * sx * sy * sy;
}

// Ghost-mirrored accessors, line-for-line the same branches as
// StaggeredCavityBase3D::uAt/vAt/wAt/nutAt. Deliberately no bounds check
// on each function's own staggered direction (i for u, j for v, k for w)
// -- the CPU originals don't have one either, and every call site in the
// Central-reduced kernels below only ever passes an in-range own-
// direction index by construction of the calling loops (verified against
// StaggeredCavityBase3D.cpp before writing this file, not assumed).
__device__ __forceinline__ double deviceUAt(const double* u, long long i, long long j, long long k,
                                             std::size_t nx, std::size_t ny, std::size_t nz, double dx,
                                             double dy, double lengthX, double lengthY, double lidVelocity) {
    const auto ii = static_cast<std::size_t>(i);
    if (j < 0) {
        return -u[deviceIndexU(ii, 0, static_cast<std::size_t>(k), nx, ny)];
    }
    if (j >= static_cast<long long>(ny)) {
        return -u[deviceIndexU(ii, ny - 1, static_cast<std::size_t>(k), nx, ny)];
    }
    if (k < 0) {
        return -u[deviceIndexU(ii, static_cast<std::size_t>(j), 0, nx, ny)];
    }
    if (k >= static_cast<long long>(nz)) {
        const double x = static_cast<double>(i) * dx;
        const double y = (static_cast<double>(j) + 0.5) * dy;
        return 2.0 * deviceLidVelocityAt(x, y, lengthX, lengthY, lidVelocity) -
               u[deviceIndexU(ii, static_cast<std::size_t>(j), nz - 1, nx, ny)];
    }
    return u[deviceIndexU(ii, static_cast<std::size_t>(j), static_cast<std::size_t>(k), nx, ny)];
}

__device__ __forceinline__ double deviceVAt(const double* v, long long i, long long j, long long k,
                                             std::size_t nx, std::size_t ny, std::size_t nz) {
    const auto jj = static_cast<std::size_t>(j);
    if (i < 0) {
        return -v[deviceIndexV(0, jj, static_cast<std::size_t>(k), nx, ny)];
    }
    if (i >= static_cast<long long>(nx)) {
        return -v[deviceIndexV(nx - 1, jj, static_cast<std::size_t>(k), nx, ny)];
    }
    if (k < 0) {
        return -v[deviceIndexV(static_cast<std::size_t>(i), jj, 0, nx, ny)];
    }
    if (k >= static_cast<long long>(nz)) {
        return -v[deviceIndexV(static_cast<std::size_t>(i), jj, nz - 1, nx, ny)];
    }
    return v[deviceIndexV(static_cast<std::size_t>(i), jj, static_cast<std::size_t>(k), nx, ny)];
}

__device__ __forceinline__ double deviceWAt(const double* w, long long i, long long j, long long k,
                                             std::size_t nx, std::size_t ny, std::size_t nz) {
    const auto kk = static_cast<std::size_t>(k);
    if (i < 0) {
        return -w[deviceIndexW(0, static_cast<std::size_t>(j), kk, nx, ny)];
    }
    if (i >= static_cast<long long>(nx)) {
        return -w[deviceIndexW(nx - 1, static_cast<std::size_t>(j), kk, nx, ny)];
    }
    if (j < 0) {
        return -w[deviceIndexW(static_cast<std::size_t>(i), 0, kk, nx, ny)];
    }
    if (j >= static_cast<long long>(ny)) {
        return -w[deviceIndexW(static_cast<std::size_t>(i), ny - 1, kk, nx, ny)];
    }
    return w[deviceIndexW(static_cast<std::size_t>(i), static_cast<std::size_t>(j), kk, nx, ny)];
}

__device__ __forceinline__ double deviceNutAt(const double* nut, long long i, long long j, long long k,
                                               std::size_t nx, std::size_t ny, std::size_t nz) {
    if (i < 0 || i >= static_cast<long long>(nx) || j < 0 || j >= static_cast<long long>(ny) || k < 0 ||
        k >= static_cast<long long>(nz)) {
        return 0.0;
    }
    return nut[deviceIndexP(static_cast<std::size_t>(i), static_cast<std::size_t>(j), static_cast<std::size_t>(k),
                             nx, ny)];
}

// Central-scheme-reduced u-momentum: StaggeredCavityBase3D::computeMomentumPredictor's
// u-loop (engine/solver/src/StaggeredCavityBase3D.cpp:176-239) with every
// schemeTransportValue(...) call replaced by its own first argument --
// exactly what that function evaluates to for ConvectionScheme::Central,
// its own body's very first branch. One thread per u-face; boundary faces
// (i==0 or i==nx) return immediately, preserved by predict()'s own
// device-to-device copy before this kernel ever runs.
__global__ void uMomentumKernel(const double* u, const double* v, const double* w, const double* nut,
                                 double* uStar, std::size_t nx, std::size_t ny, std::size_t nz, double dx,
                                 double dy, double dz, double lengthX, double lengthY, double viscosity,
                                 double lidVelocity, double dt) {
    const auto i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto j = static_cast<std::size_t>(blockIdx.y) * blockDim.y + threadIdx.y;
    const auto k = static_cast<std::size_t>(blockIdx.z) * blockDim.z + threadIdx.z;
    if (i > nx || j >= ny || k >= nz) {
        return;
    }
    if (i == 0 || i == nx) {
        return;
    }
    const auto li = static_cast<long long>(i);
    const auto lj = static_cast<long long>(j);
    const auto lk = static_cast<long long>(k);

    const double uHere = deviceUAt(u, li, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
    const double uIp1 = deviceUAt(u, li + 1, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
    const double uIm1 = deviceUAt(u, li - 1, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
    const double uCenterI = 0.5 * (uHere + uIp1);
    const double uCenterIm1 = 0.5 * (uIm1 + uHere);
    const double duudx = (uCenterI * uCenterI - uCenterIm1 * uCenterIm1) / dx;

    const double uJp1 = deviceUAt(u, li, lj + 1, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
    const double uJm1 = deviceUAt(u, li, lj - 1, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
    const double uEdgeJp = 0.5 * (uHere + uJp1);
    const double uEdgeJm = 0.5 * (uJm1 + uHere);
    const double vEdgeIJp =
        0.5 * (deviceVAt(v, li - 1, lj + 1, lk, nx, ny, nz) + deviceVAt(v, li, lj + 1, lk, nx, ny, nz));
    const double vEdgeIJ = 0.5 * (deviceVAt(v, li - 1, lj, lk, nx, ny, nz) + deviceVAt(v, li, lj, lk, nx, ny, nz));
    const double duvdy = (vEdgeIJp * uEdgeJp - vEdgeIJ * uEdgeJm) / dy;

    const double uKp1 = deviceUAt(u, li, lj, lk + 1, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
    const double uKm1 = deviceUAt(u, li, lj, lk - 1, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity);
    const double uEdgeKp = 0.5 * (uHere + uKp1);
    const double uEdgeKm = 0.5 * (uKm1 + uHere);
    const double wEdgeIKp =
        0.5 * (deviceWAt(w, li - 1, lj, lk + 1, nx, ny, nz) + deviceWAt(w, li, lj, lk + 1, nx, ny, nz));
    const double wEdgeIK = 0.5 * (deviceWAt(w, li - 1, lj, lk, nx, ny, nz) + deviceWAt(w, li, lj, lk, nx, ny, nz));
    const double duwdz = (wEdgeIKp * uEdgeKp - wEdgeIK * uEdgeKm) / dz;

    const double gammaE = viscosity + deviceNutAt(nut, li, lj, lk, nx, ny, nz);
    const double gammaW = viscosity + deviceNutAt(nut, li - 1, lj, lk, nx, ny, nz);
    const double gammaTransverse =
        viscosity + 0.5 * (deviceNutAt(nut, li - 1, lj, lk, nx, ny, nz) + deviceNutAt(nut, li, lj, lk, nx, ny, nz));

    const double diffusionU = (gammaE * (uIp1 - uHere) - gammaW * (uHere - uIm1)) / (dx * dx) +
                               gammaTransverse * (uJp1 - 2.0 * uHere + uJm1) / (dy * dy) +
                               gammaTransverse * (uKp1 - 2.0 * uHere + uKm1) / (dz * dz);

    uStar[deviceIndexU(i, j, k, nx, ny)] = uHere + dt * (-(duudx + duvdy + duwdz) + diffusionU);
}

// Central-scheme-reduced v-momentum: StaggeredCavityBase3D.cpp:241-301.
// Own direction j; cross terms with w (z) and u (x).
__global__ void vMomentumKernel(const double* u, const double* v, const double* w, const double* nut,
                                 double* vStar, std::size_t nx, std::size_t ny, std::size_t nz, double dx,
                                 double dy, double dz, double lengthX, double lengthY, double viscosity,
                                 double lidVelocity, double dt) {
    const auto i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto j = static_cast<std::size_t>(blockIdx.y) * blockDim.y + threadIdx.y;
    const auto k = static_cast<std::size_t>(blockIdx.z) * blockDim.z + threadIdx.z;
    if (i >= nx || j > ny || k >= nz) {
        return;
    }
    if (j == 0 || j == ny) {
        return;
    }
    const auto li = static_cast<long long>(i);
    const auto lj = static_cast<long long>(j);
    const auto lk = static_cast<long long>(k);

    const double vHere = deviceVAt(v, li, lj, lk, nx, ny, nz);
    const double vJp1 = deviceVAt(v, li, lj + 1, lk, nx, ny, nz);
    const double vJm1 = deviceVAt(v, li, lj - 1, lk, nx, ny, nz);
    const double vCenterJ = 0.5 * (vHere + vJp1);
    const double vCenterJm1 = 0.5 * (vJm1 + vHere);
    const double dvvdy = (vCenterJ * vCenterJ - vCenterJm1 * vCenterJm1) / dy;

    const double vKp1 = deviceVAt(v, li, lj, lk + 1, nx, ny, nz);
    const double vKm1 = deviceVAt(v, li, lj, lk - 1, nx, ny, nz);
    const double vEdgeKp = 0.5 * (vHere + vKp1);
    const double vEdgeKm = 0.5 * (vKm1 + vHere);
    const double wEdgeJKp =
        0.5 * (deviceWAt(w, li, lj - 1, lk + 1, nx, ny, nz) + deviceWAt(w, li, lj, lk + 1, nx, ny, nz));
    const double wEdgeJK = 0.5 * (deviceWAt(w, li, lj - 1, lk, nx, ny, nz) + deviceWAt(w, li, lj, lk, nx, ny, nz));
    const double dvwdz = (wEdgeJKp * vEdgeKp - wEdgeJK * vEdgeKm) / dz;

    const double vIp1 = deviceVAt(v, li + 1, lj, lk, nx, ny, nz);
    const double vIm1 = deviceVAt(v, li - 1, lj, lk, nx, ny, nz);
    const double vEdgeIp = 0.5 * (vHere + vIp1);
    const double vEdgeIm = 0.5 * (vIm1 + vHere);
    const double uEdgeJIp = 0.5 * (deviceUAt(u, li + 1, lj - 1, lk, nx, ny, nz, dx, dy, lengthX, lengthY,
                                              lidVelocity) +
                                    deviceUAt(u, li + 1, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity));
    const double uEdgeJI = 0.5 * (deviceUAt(u, li, lj - 1, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity) +
                                   deviceUAt(u, li, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity));
    const double dvudx = (uEdgeJIp * vEdgeIp - uEdgeJI * vEdgeIm) / dx;

    const double gammaN = viscosity + deviceNutAt(nut, li, lj, lk, nx, ny, nz);
    const double gammaS = viscosity + deviceNutAt(nut, li, lj - 1, lk, nx, ny, nz);
    const double gammaTransverse =
        viscosity + 0.5 * (deviceNutAt(nut, li, lj - 1, lk, nx, ny, nz) + deviceNutAt(nut, li, lj, lk, nx, ny, nz));

    const double diffusionV = (gammaN * (vJp1 - vHere) - gammaS * (vHere - vJm1)) / (dy * dy) +
                               gammaTransverse * (vIp1 - 2.0 * vHere + vIm1) / (dx * dx) +
                               gammaTransverse * (vKp1 - 2.0 * vHere + vKm1) / (dz * dz);

    vStar[deviceIndexV(i, j, k, nx, ny)] = vHere + dt * (-(dvvdy + dvwdz + dvudx) + diffusionV);
}

// Central-scheme-reduced w-momentum: StaggeredCavityBase3D.cpp:303-363.
// Own direction k; cross terms with u (x) and v (y).
__global__ void wMomentumKernel(const double* u, const double* v, const double* w, const double* nut,
                                 double* wStar, std::size_t nx, std::size_t ny, std::size_t nz, double dx,
                                 double dy, double dz, double lengthX, double lengthY, double viscosity,
                                 double lidVelocity, double dt) {
    const auto i = static_cast<std::size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    const auto j = static_cast<std::size_t>(blockIdx.y) * blockDim.y + threadIdx.y;
    const auto k = static_cast<std::size_t>(blockIdx.z) * blockDim.z + threadIdx.z;
    if (i >= nx || j >= ny || k > nz) {
        return;
    }
    if (k == 0 || k == nz) {
        return;
    }
    const auto li = static_cast<long long>(i);
    const auto lj = static_cast<long long>(j);
    const auto lk = static_cast<long long>(k);

    const double wHere = deviceWAt(w, li, lj, lk, nx, ny, nz);
    const double wKp1 = deviceWAt(w, li, lj, lk + 1, nx, ny, nz);
    const double wKm1 = deviceWAt(w, li, lj, lk - 1, nx, ny, nz);
    const double wCenterK = 0.5 * (wHere + wKp1);
    const double wCenterKm1 = 0.5 * (wKm1 + wHere);
    const double dwwdz = (wCenterK * wCenterK - wCenterKm1 * wCenterKm1) / dz;

    const double wIp1 = deviceWAt(w, li + 1, lj, lk, nx, ny, nz);
    const double wIm1 = deviceWAt(w, li - 1, lj, lk, nx, ny, nz);
    const double wEdgeIp = 0.5 * (wHere + wIp1);
    const double wEdgeIm = 0.5 * (wIm1 + wHere);
    const double uEdgeKIp = 0.5 * (deviceUAt(u, li + 1, lj, lk - 1, nx, ny, nz, dx, dy, lengthX, lengthY,
                                              lidVelocity) +
                                    deviceUAt(u, li + 1, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity));
    const double uEdgeKI = 0.5 * (deviceUAt(u, li, lj, lk - 1, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity) +
                                   deviceUAt(u, li, lj, lk, nx, ny, nz, dx, dy, lengthX, lengthY, lidVelocity));
    const double dwudx = (uEdgeKIp * wEdgeIp - uEdgeKI * wEdgeIm) / dx;

    const double wJp1 = deviceWAt(w, li, lj + 1, lk, nx, ny, nz);
    const double wJm1 = deviceWAt(w, li, lj - 1, lk, nx, ny, nz);
    const double wEdgeJp = 0.5 * (wHere + wJp1);
    const double wEdgeJm = 0.5 * (wJm1 + wHere);
    const double vEdgeKJp =
        0.5 * (deviceVAt(v, li, lj + 1, lk - 1, nx, ny, nz) + deviceVAt(v, li, lj + 1, lk, nx, ny, nz));
    const double vEdgeKJ = 0.5 * (deviceVAt(v, li, lj, lk - 1, nx, ny, nz) + deviceVAt(v, li, lj, lk, nx, ny, nz));
    const double dwvdy = (vEdgeKJp * wEdgeJp - vEdgeKJ * wEdgeJm) / dy;

    const double gammaF = viscosity + deviceNutAt(nut, li, lj, lk, nx, ny, nz);
    const double gammaB = viscosity + deviceNutAt(nut, li, lj, lk - 1, nx, ny, nz);
    const double gammaTransverse =
        viscosity + 0.5 * (deviceNutAt(nut, li, lj, lk - 1, nx, ny, nz) + deviceNutAt(nut, li, lj, lk, nx, ny, nz));

    const double diffusionW = (gammaF * (wKp1 - wHere) - gammaB * (wHere - wKm1)) / (dz * dz) +
                               gammaTransverse * (wIp1 - 2.0 * wHere + wIm1) / (dx * dx) +
                               gammaTransverse * (wJp1 - 2.0 * wHere + wJm1) / (dy * dy);

    wStar[deviceIndexW(i, j, k, nx, ny)] = wHere + dt * (-(dwwdz + dwudx + dwvdy) + diffusionW);
}

} // namespace

struct MomentumPredictorCuda::Impl {
    double* u = nullptr;
    double* v = nullptr;
    double* w = nullptr;
    double* nut = nullptr;
    double* uStar = nullptr;
    double* vStar = nullptr;
    double* wStar = nullptr;

    ~Impl() {
        cudaFree(u);
        cudaFree(v);
        cudaFree(w);
        cudaFree(nut);
        cudaFree(uStar);
        cudaFree(vStar);
        cudaFree(wStar);
    }
};

MomentumPredictorCuda::MomentumPredictorCuda(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                              double lengthY, double lengthZ, double viscosity,
                                              double lidVelocity)
    : impl_(std::make_unique<Impl>()), nx_(nx), ny_(ny), nz_(nz), lengthX_(lengthX), lengthY_(lengthY),
      lengthZ_(lengthZ), dx_(lengthX / static_cast<double>(nx)), dy_(lengthY / static_cast<double>(ny)),
      dz_(lengthZ / static_cast<double>(nz)), viscosity_(viscosity), lidVelocity_(lidVelocity),
      available_(false) {
    int deviceCount = 0;
    if (cudaGetDeviceCount(&deviceCount) != cudaSuccess || deviceCount == 0) {
        return; // available_ stays false; impl_'s still-null buffers are freed harmlessly by ~Impl()
    }

    const std::size_t uSize = (nx_ + 1) * ny_ * nz_;
    const std::size_t vSize = nx_ * (ny_ + 1) * nz_;
    const std::size_t wSize = nx_ * ny_ * (nz_ + 1);
    const std::size_t pSize = nx_ * ny_ * nz_;

    // Every allocation checked individually -- seven persistent buffers,
    // same reasoning as ConjugateGradientSolverCuda's own constructor: a
    // large-grid VRAM exhaustion is a real failure mode, and impl_ is
    // already a fully-constructed unique_ptr<Impl> at this point, so a
    // throw here still unwinds through ~Impl(), freeing whatever subset
    // already succeeded.
    const auto checkedMalloc = [](double** ptr, std::size_t count) {
        const cudaError_t status = cudaMalloc(ptr, count * sizeof(double));
        if (status != cudaSuccess) {
            throw std::runtime_error(std::string("MomentumPredictorCuda: cudaMalloc failed: ") +
                                      cudaGetErrorString(status));
        }
    };
    checkedMalloc(&impl_->u, uSize);
    checkedMalloc(&impl_->v, vSize);
    checkedMalloc(&impl_->w, wSize);
    checkedMalloc(&impl_->nut, pSize);
    checkedMalloc(&impl_->uStar, uSize);
    checkedMalloc(&impl_->vStar, vSize);
    checkedMalloc(&impl_->wStar, wSize);

    available_ = true;
}

MomentumPredictorCuda::~MomentumPredictorCuda() = default;

MomentumPredictorCuda::Result MomentumPredictorCuda::predict(const std::vector<double>& u,
                                                               const std::vector<double>& v,
                                                               const std::vector<double>& w,
                                                               const std::vector<double>& nut, double dt) const {
    if (!available_) {
        throw std::runtime_error("MomentumPredictorCuda::predict: no CUDA-capable device available");
    }
    const std::size_t uSize = (nx_ + 1) * ny_ * nz_;
    const std::size_t vSize = nx_ * (ny_ + 1) * nz_;
    const std::size_t wSize = nx_ * ny_ * (nz_ + 1);
    const std::size_t pSize = nx_ * ny_ * nz_;
    if (u.size() != uSize || v.size() != vSize || w.size() != wSize) {
        throw std::invalid_argument("MomentumPredictorCuda::predict: u/v/w size mismatch");
    }
    if (!nut.empty() && nut.size() != pSize) {
        throw std::invalid_argument("MomentumPredictorCuda::predict: nut size mismatch");
    }

    Impl& impl = *impl_;
    cudaMemcpy(impl.u, u.data(), uSize * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(impl.v, v.data(), vSize * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(impl.w, w.data(), wSize * sizeof(double), cudaMemcpyHostToDevice);
    // nut is reused across calls -- a laminar call must explicitly
    // re-zero it, since a prior turbulent call could otherwise leave
    // stale values behind (it cannot rely on "it was zero at
    // construction" once any later call has dirtied it).
    if (nut.empty()) {
        cudaMemset(impl.nut, 0, pSize * sizeof(double));
    } else {
        cudaMemcpy(impl.nut, nut.data(), pSize * sizeof(double), cudaMemcpyHostToDevice);
    }

    // Boundary faces preserved by construction, not by kernel special-case:
    // copy the whole field first, then the kernels below only ever
    // overwrite interior faces -- matching the CPU contract exactly
    // (callers pass copies of u_/v_/w_ as uStar/vStar/wStar, and
    // computeMomentumPredictor only ever writes interior faces).
    cudaMemcpy(impl.uStar, impl.u, uSize * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(impl.vStar, impl.v, vSize * sizeof(double), cudaMemcpyDeviceToDevice);
    cudaMemcpy(impl.wStar, impl.w, wSize * sizeof(double), cudaMemcpyDeviceToDevice);

    const dim3 blockDim3(8, 8, 8);
    const dim3 gridDimU(static_cast<unsigned>((nx_ + 1 + 7) / 8), static_cast<unsigned>((ny_ + 7) / 8),
                         static_cast<unsigned>((nz_ + 7) / 8));
    const dim3 gridDimV(static_cast<unsigned>((nx_ + 7) / 8), static_cast<unsigned>((ny_ + 1 + 7) / 8),
                         static_cast<unsigned>((nz_ + 7) / 8));
    const dim3 gridDimW(static_cast<unsigned>((nx_ + 7) / 8), static_cast<unsigned>((ny_ + 7) / 8),
                         static_cast<unsigned>((nz_ + 1 + 7) / 8));

    // No dependency between the three launches -- each writes a disjoint
    // output buffer and only reads from the already fully-uploaded
    // u/v/w/nut inputs -- so no stream synchronization is needed between
    // them, only before the final device->host copy below.
    uMomentumKernel<<<gridDimU, blockDim3>>>(impl.u, impl.v, impl.w, impl.nut, impl.uStar, nx_, ny_, nz_, dx_,
                                              dy_, dz_, lengthX_, lengthY_, viscosity_, lidVelocity_, dt);
    vMomentumKernel<<<gridDimV, blockDim3>>>(impl.u, impl.v, impl.w, impl.nut, impl.vStar, nx_, ny_, nz_, dx_,
                                              dy_, dz_, lengthX_, lengthY_, viscosity_, lidVelocity_, dt);
    wMomentumKernel<<<gridDimW, blockDim3>>>(impl.u, impl.v, impl.w, impl.nut, impl.wStar, nx_, ny_, nz_, dx_,
                                              dy_, dz_, lengthX_, lengthY_, viscosity_, lidVelocity_, dt);
    cudaDeviceSynchronize();

    Result result;
    result.uStar.resize(uSize);
    result.vStar.resize(vSize);
    result.wStar.resize(wSize);
    cudaMemcpy(result.uStar.data(), impl.uStar, uSize * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(result.vStar.data(), impl.vStar, vSize * sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(result.wStar.data(), impl.wStar, wSize * sizeof(double), cudaMemcpyDeviceToHost);

    const cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        throw std::runtime_error(std::string("MomentumPredictorCuda::predict: CUDA error: ") +
                                  cudaGetErrorString(status));
    }
    return result;
}

} // namespace aether::gpu
