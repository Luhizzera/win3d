#include "aether/solver/StaggeredLidDrivenCavitySolver3D.hpp"

namespace aether::solver {

StaggeredLidDrivenCavitySolver3D::StaggeredLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX, double lengthY,
             double lengthZ, double viscosity, double lidVelocity, ConvectionScheme convection)
    : StaggeredCavityBase3D(nx, ny, nz, lengthX, lengthY, lengthZ, viscosity, lidVelocity, convection) {}

void StaggeredLidDrivenCavitySolver3D::step(double dt) {
    std::vector<double> uStar = u_;
    std::vector<double> vStar = v_;
    std::vector<double> wStar = w_;

    computeMomentumPredictor(uStar, vStar, wStar, dt);
    projectToDivergenceFree(uStar, vStar, wStar, dt);
    time_ += dt;
}

} // namespace aether::solver
