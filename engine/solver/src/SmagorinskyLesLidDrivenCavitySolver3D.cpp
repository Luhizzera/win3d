#include "aether/solver/SmagorinskyLesLidDrivenCavitySolver3D.hpp"

#include <algorithm>
#include <cmath>

namespace aether::solver {

namespace {
// The same von Karman constant every other closure in this project uses;
// here it appears only in the near-wall cap on the mixing length (see the
// class header on why that stands in for van Driest damping).
constexpr double kKarman = 0.41;
} // namespace

SmagorinskyLesLidDrivenCavitySolver3D::SmagorinskyLesLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX, double lengthY,
             double lengthZ, double viscosity, double lidVelocity, double smagorinskyConstant)
    : StaggeredCavityBase3D(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity),
      smagorinskyConstant_(smagorinskyConstant), nut_(nx * ny * nz, 0.0) {
    setEddyViscosityField(&nut_);
}

double SmagorinskyLesLidDrivenCavitySolver3D::subgridViscosity(std::size_t i, std::size_t j, std::size_t k) const {
    return nut_[indexP(i, j, k)];
}

double SmagorinskyLesLidDrivenCavitySolver3D::filterWidth() const {
    return std::cbrt(dx_ * dy_ * dz_);
}

void SmagorinskyLesLidDrivenCavitySolver3D::updateSubgridViscosity() {
    const double delta = filterWidth();

    for (std::size_t k = 0; k < nz_; ++k) {
        for (std::size_t j = 0; j < ny_; ++j) {
            for (std::size_t i = 0; i < nx_; ++i) {
                const auto li = static_cast<long long>(i);
                const auto lj = static_cast<long long>(j);
                const auto lk = static_cast<long long>(k);

                // Identical strain-rate construction to
                // MixingLengthLidDrivenCavitySolver3D: own-axis components
                // are exact face differences; cross components average the
                // two central differences taken at the cell's two flanking
                // faces along the component's own axis. Reused rather than
                // reinvented -- see that class's header for the derivation.
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

                // l_s = min(Cs*Delta, kappa*d_wall): the Smagorinsky
                // length is set by the MESH (Cs*Delta), capped near a wall
                // by the geometric mixing length so nu_sgs vanishes there
                // -- the deliberate stand-in for van Driest damping
                // explained in the class header.
                const double lengthScale =
                    std::min(smagorinskyConstant_ * delta, kKarman * wallDistanceAt(i, j, k));

                nut_[indexP(i, j, k)] = lengthScale * lengthScale * strainMagnitude;
            }
        }
    }
}

void SmagorinskyLesLidDrivenCavitySolver3D::step(double dt) {
    updateSubgridViscosity();

    std::vector<double> uStar = u_;
    std::vector<double> vStar = v_;
    std::vector<double> wStar = w_;

    computeMomentumPredictor(uStar, vStar, wStar, dt);
    projectToDivergenceFree(uStar, vStar, wStar, dt);
    time_ += dt;
}

} // namespace aether::solver
