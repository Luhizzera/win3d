#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// The first genuinely 3D turbulence-convection coupling in this project:
// StaggeredLidDrivenCavitySolver3D's solid-wall MAC-grid box, closed with
// the same Prandtl mixing-length eddy-viscosity closure already validated
// on the 2D cavity (MixingLengthLidDrivenCavitySolver2D). Every turbulence
// closure so far that saw a real 2D/3D convecting field
// (MixingLengthLidDrivenCavitySolver2D, KEpsilonLidDrivenCavitySolver2D,
// KOmegaSSTLidDrivenCavitySolver2D) was still confined to a flat, colocated
// 2D grid; this is the first one on a staggered 3D grid.
//
// nu_t(x,y,z) = l_m(x,y,z)^2 * |S(x,y,z)|, with:
//   l_m = min(kappa * d_wall, 0.09 * min(lengthX,lengthY,lengthZ)/2) -- the
//     same von-Karman-near-wall / Escudier-cap formula as every other
//     mixing-length closure in this project, now with d_wall = distance to
//     the nearest of all *six* solid walls.
//   |S| = sqrt(2*Sij*Sij) = sqrt(2*(Sxx^2+Syy^2+Szz^2) + 4*(Sxy^2+Sxz^2+Syz^2)),
//     the full 3D strain-rate magnitude (the 2D cavity's formula is the
//     special case Szz=Sxz=Syz=0). Computed at cell centers (the same
//     points pressure lives at) by combining a strain component's *own*
//     axis (an exact face-difference: e.g. dudx from u's two flanking
//     x-faces of the cell) with its *cross* components (e.g. dudy: averaged
//     from the two du/dy values evaluated at the cell's two x-faces, each a
//     central difference in y across that face's own ghost-mirrored
//     neighbors) -- the same edge-averaging pattern step()'s own convective
//     cross-terms (duvdy, dvudx, ...) already use on this staggered grid,
//     reused here instead of inventing a different interpolation scheme.
//
// nu_t is folded into each component's momentum diffusion as a *locally
// sampled* effective viscosity: exact face-varying (nu + nu_t at the two
// flanking pressure cells) along the component's own axis -- e.g. u's
// x-direction term uses nu_t at cells (i-1,j,k) and (i,j,k), the two cells
// u(i,j,k) genuinely sits between -- and a single cell-averaged value along
// the two transverse axes (e.g. u's y- and z-direction terms both use
// 0.5*(nut(i-1,j,k)+nut(i,j,k)), the same pair, rather than resolving a
// separate value per transverse sub-face). This is a further simplification
// of the already-simplified "quasi-laminar form" MixingLengthLidDrivenCavitySolver2D
// documents (which itself omits the full variable-viscosity stress-tensor
// cross terms) -- justified the same way: a fully rigorous treatment would
// need nu_t resolved at every edge of the staggered grid (12 distinct
// locations per cell in 3D, each needing its own 4-cell average), which is
// real additional bookkeeping for an accuracy gain not claimed or needed by
// this project's validation philosophy (structural/first-principles checks,
// not literature benchmark matching). nu_t is exactly zero at all six walls
// by construction (mixing length vanishes there, and nutAt() returns 0.0 for
// any out-of-domain cell).
//
// Same ghost-mirror walls, tapered lid, matrix-free CG pressure projection,
// and "no literature benchmark" validation as StaggeredLidDrivenCavitySolver3D:
// divergence stays bounded, the primary-vortex topology still holds, and
// nu_t is checked to vanish at the walls and be non-negative.
class MixingLengthLidDrivenCavitySolver3D {
public:
    MixingLengthLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                         double lengthY, double lengthZ, double viscosity,
                                         double lidVelocity);

    // Conservative estimate from molecular viscosity and the current
    // (possibly still-zero, if called before any step()) eddy viscosity --
    // same caveat as MixingLengthLidDrivenCavitySolver2D::stableTimeStep().
    double stableTimeStep() const;

    void step(double dt);

    double u(std::size_t i, std::size_t j, std::size_t k) const; // i in [0,nx]
    double v(std::size_t i, std::size_t j, std::size_t k) const; // j in [0,ny]
    double w(std::size_t i, std::size_t j, std::size_t k) const; // k in [0,nz]
    double pressure(std::size_t i, std::size_t j, std::size_t k) const;
    double eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const;
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
    double wallDistanceAt(std::size_t i, std::size_t j, std::size_t k) const;

    double uAt(long long i, long long j, long long k) const;
    double vAt(long long i, long long j, long long k) const;
    double wAt(long long i, long long j, long long k) const;
    double pAt(long long i, long long j, long long k) const;
    // Cell-centered eddy viscosity; 0.0 for any cell index outside
    // [0,nx)x[0,ny)x[0,nz) (nu_t vanishes at/beyond a solid wall).
    double nutAt(long long i, long long j, long long k) const;

    void updateEddyViscosity();

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
    std::vector<double> nut_;
    double time_ = 0.0;
};

} // namespace aether::solver
