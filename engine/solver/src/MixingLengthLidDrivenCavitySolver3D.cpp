#include "aether/solver/MixingLengthLidDrivenCavitySolver3D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

namespace {
constexpr double kKarman = 0.41;
constexpr double kEscudierFactor = 0.09;
} // namespace

MixingLengthLidDrivenCavitySolver3D::MixingLengthLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX, double lengthY,
             double lengthZ, double viscosity, double lidVelocity, bool useGpu)
    : StaggeredCavityBase3D(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity,
                            ConvectionScheme::Central, useGpu),
      nut_(nx * ny * nz, 0.0) {
    // Lets the base's momentum predictor and time-step limit see this
    // closure's nu_t without a virtual call in the inner loop.
    setEddyViscosityField(&nut_);
}

double MixingLengthLidDrivenCavitySolver3D::eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const {
    return nut_[indexP(i, j, k)];
}

void MixingLengthLidDrivenCavitySolver3D::updateEddyViscosity() {
    const double lengthCap = kEscudierFactor * std::min({lengthX_, lengthY_, lengthZ_}) / 2.0;

    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                const double dudx = (uAt(li + 1, lj, lk) - uAt(li, lj, lk)) / dx_;
                const double dvdy = (vAt(li, lj + 1, lk) - vAt(li, lj, lk)) / dy_;
                const double dwdz = (wAt(li, lj, lk + 1) - wAt(li, lj, lk)) / dz_;

                const double dudyFaceI = (uAt(li, lj + 1, lk) - uAt(li, lj - 1, lk)) / (2.0 * dy_);
                const double dudyFaceIp1 = (uAt(li + 1, lj + 1, lk) - uAt(li + 1, lj - 1, lk)) / (2.0 * dy_);
                const double dudy = 0.5 * (dudyFaceI + dudyFaceIp1);

                const double dvdxFaceJ = (vAt(li + 1, lj, lk) - vAt(li - 1, lj, lk)) / (2.0 * dx_);
                const double dvdxFaceJp1 = (vAt(li + 1, lj + 1, lk) - vAt(li - 1, lj + 1, lk)) / (2.0 * dx_);
                const double dvdx = 0.5 * (dvdxFaceJ + dvdxFaceJp1);

                const double dudzFaceI = (uAt(li, lj, lk + 1) - uAt(li, lj, lk - 1)) / (2.0 * dz_);
                const double dudzFaceIp1 = (uAt(li + 1, lj, lk + 1) - uAt(li + 1, lj, lk - 1)) / (2.0 * dz_);
                const double dudz = 0.5 * (dudzFaceI + dudzFaceIp1);

                const double dwdxFaceK = (wAt(li + 1, lj, lk) - wAt(li - 1, lj, lk)) / (2.0 * dx_);
                const double dwdxFaceKp1 = (wAt(li + 1, lj, lk + 1) - wAt(li - 1, lj, lk + 1)) / (2.0 * dx_);
                const double dwdx = 0.5 * (dwdxFaceK + dwdxFaceKp1);

                const double dvdzFaceJ = (vAt(li, lj, lk + 1) - vAt(li, lj, lk - 1)) / (2.0 * dz_);
                const double dvdzFaceJp1 = (vAt(li, lj + 1, lk + 1) - vAt(li, lj + 1, lk - 1)) / (2.0 * dz_);
                const double dvdz = 0.5 * (dvdzFaceJ + dvdzFaceJp1);

                const double dwdyFaceK = (wAt(li, lj + 1, lk) - wAt(li, lj - 1, lk)) / (2.0 * dy_);
                const double dwdyFaceKp1 = (wAt(li, lj + 1, lk + 1) - wAt(li, lj - 1, lk + 1)) / (2.0 * dy_);
                const double dwdy = 0.5 * (dwdyFaceK + dwdyFaceKp1);

                const double sxy = 0.5 * (dudy + dvdx);
                const double sxz = 0.5 * (dudz + dwdx);
                const double syz = 0.5 * (dvdz + dwdy);

                const double strainMagnitude = std::sqrt(2.0 * (dudx * dudx + dvdy * dvdy + dwdz * dwdz) +
                                                           4.0 * (sxy * sxy + sxz * sxz + syz * syz));

                const double wallDist = wallDistanceAt(i, j, k);
                const double mixingLength = std::min(kKarman * wallDist, lengthCap);

                nut_[indexP(i, j, k)] = mixingLength * mixingLength * strainMagnitude;
            }
        }
    }
}

void MixingLengthLidDrivenCavitySolver3D::step(double dt) {
    updateEddyViscosity();

    std::vector<double> uStar = u_;
    std::vector<double> vStar = v_;
    std::vector<double> wStar = w_;

    computeMomentumPredictor(uStar, vStar, wStar, dt);
    projectToDivergenceFree(uStar, vStar, wStar, dt);
    time_ += dt;
}

} // namespace aether::solver
