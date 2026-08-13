#pragma once

#include "aether/solver/StaggeredCavityBase3D.hpp"

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
// The staggered-grid plumbing (indexing, ghost mirrors, momentum predictor,
// pressure projection, stableTimeStep, maxDivergence) lives in
// StaggeredCavityBase3D; only the closure itself is here.
class MixingLengthLidDrivenCavitySolver3D : public StaggeredCavityBase3D {
public:
    MixingLengthLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                         double lengthY, double lengthZ, double viscosity,
                                         double lidVelocity);

    void step(double dt);

    double eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const;

private:
    void updateEddyViscosity();

    std::vector<double> nut_;
};

} // namespace aether::solver
