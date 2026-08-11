#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Solves the 3D incompressible Navier-Stokes equations on a box with solid
// walls on all six faces and a moving lid on top (z = lengthZ) -- the 3D
// generalization of LidDrivenCavitySolver2D, on the staggered (MAC) grid
// StaggeredNavierStokesSolver3D introduced for the periodic case. This
// closes the "solid-wall 3D" gap flagged when that class shipped.
//
// Unlike the periodic solver (where u, v, w each have nx*ny*nz entries,
// since a periodic face coincides with cell 0), a solid-walled domain gives
// each velocity component genuine boundary FACES that are not the same as
// any interior face, so the arrays are sized one larger along their own
// axis: u has (nx+1)*ny*nz entries (i in [0,nx]), v has nx*(ny+1)*nz
// (j in [0,ny]), w has nx*ny*(nz+1) (k in [0,nz]). The two array entries at
// each end of a component's own axis (e.g. u(0,j,k) and u(nx,j,k)) are the
// actual physical boundary values themselves -- direct Dirichlet fixes
// (always 0: no flow penetrates any solid wall, and the lid moves only
// tangentially), never solved for and never touched after construction.
//
// Every *other* direction is a component's tangential direction relative
// to a wall it doesn't have a face on (e.g. u has no face at y=0/lengthY
// or z=0/lengthZ) -- there, boundary conditions are enforced by Dirichlet
// ghost-mirroring (ghost = 2*wallValue - interior), the same technique
// LidDrivenCavitySolver2D and MixingLengthLidDrivenCavitySolver2D already
// use. The lid's tapered tangential velocity
// (lidVelocity * sin^2(pi*x/lengthX) * sin^2(pi*y/lengthY) -- tapering to
// zero along all four edges of the top face, the 3D generalization of the
// 2D solver's corner-only taper, for the same reason: an actually
// discontinuous lid is a genuine pressure singularity) only ever enters
// through u's ghost mirror at the top wall (k >= nz); every other ghost
// mirror in this class is homogeneous (wall value 0), since only the u
// component is tangential to the lid's own direction of motion.
//
// Same Chorin projection (matrix-free Conjugate Gradient pressure solve,
// one cell pinned to remove the Neumann null space) and same "no
// literature benchmark" validation philosophy as every other Navier-Stokes
// solver in this project.
class StaggeredLidDrivenCavitySolver3D {
public:
    StaggeredLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                      double lengthY, double lengthZ, double viscosity,
                                      double lidVelocity);

    double stableTimeStep() const;
    void step(double dt);

    double u(std::size_t i, std::size_t j, std::size_t k) const; // i in [0,nx]
    double v(std::size_t i, std::size_t j, std::size_t k) const; // j in [0,ny]
    double w(std::size_t i, std::size_t j, std::size_t k) const; // k in [0,nz]
    double pressure(std::size_t i, std::size_t j, std::size_t k) const;
    double time() const { return time_; }

    double maxDivergence() const;

private:
    std::size_t indexU(std::size_t i, std::size_t j, std::size_t k) const {
        return i + j * (nx_ + 1) + k * (nx_ + 1) * ny_;
    }
    std::size_t indexV(std::size_t i, std::size_t j, std::size_t k) const {
        return i + j * nx_ + k * nx_ * (ny_ + 1);
    }
    std::size_t indexW(std::size_t i, std::size_t j, std::size_t k) const { return i + j * nx_ + k * nx_ * ny_; }
    std::size_t indexP(std::size_t i, std::size_t j, std::size_t k) const { return i + j * nx_ + k * nx_ * ny_; }

    double lidVelocityAt(double x, double y) const;

    // Direct array access for i in [0,nx]; ghost-mirrors (homogeneous,
    // except the lid at k>=nz) for j/k out of range.
    double uAt(long long i, long long j, long long k) const;
    double vAt(long long i, long long j, long long k) const;
    double wAt(long long i, long long j, long long k) const;
    // Neumann mirror (ghost = interior) across all six walls.
    double pAt(long long i, long long j, long long k) const;

    std::vector<double> applyLaplacian(const std::vector<double>& x) const;
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    void projectToDivergenceFree(std::vector<double>& uStar, std::vector<double>& vStar,
                                  std::vector<double>& wStar, double dt);

    std::size_t nx_;
    std::size_t ny_;
    std::size_t nz_;
    double lengthX_;
    double lengthY_;
    double lengthZ_;
    double dx_;
    double dy_;
    double dz_;
    double viscosity_;
    double lidVelocity_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> w_;
    std::vector<double> p_;
    double time_ = 0.0;
};

} // namespace aether::solver
