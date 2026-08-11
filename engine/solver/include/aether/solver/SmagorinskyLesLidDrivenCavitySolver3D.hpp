#pragma once

#include "aether/solver/StaggeredCavityBase3D.hpp"

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 6, next stage after the RANS closures: Large Eddy Simulation with
// the Smagorinsky subgrid-scale model, on the same solid-wall 3D staggered
// (MAC) cavity every 3D closure in this project uses
// (StaggeredLidDrivenCavitySolver3D).
//
// **The conceptual break from every earlier closure here is what the
// length scale is tied to**, and it is worth stating plainly because the
// formula otherwise looks deceptively similar to
// MixingLengthLidDrivenCavitySolver3D's:
//
//   RANS (mixing length, k-epsilon, k-omega SST): models *all* turbulent
//     motion. Its length scale comes from the FLOW GEOMETRY (wall
//     distance, or a transported turbulence scale). Refining the mesh
//     resolves the same modelled physics more accurately, but does not
//     systematically shrink nu_t -- the model is still standing in for
//     the entire turbulent spectrum.
//
//   LES (this class): resolves the large, energy-carrying eddies directly
//     and models only what falls below the grid filter width. Its length
//     scale is the FILTER WIDTH Delta, tied to the MESH:
//     Delta = cbrt(dx*dy*dz), the standard choice for an anisotropic cell.
//     So refining the mesh shrinks Delta, shrinks nu_sgs, and moves the
//     simulation toward DNS in the limit -- nu_sgs -> 0 as Delta -> 0.
//     That limit is the defining property of an LES model, and unlike
//     most claims about turbulence models it is directly measurable, so
//     this class's test suite checks it explicitly rather than asserting
//     it in a comment.
//
// nu_sgs = l_s^2 * |S|, with |S| = sqrt(2*Sij*Sij) -- the same full 3D
// strain-rate magnitude, computed by the same staggered-grid edge-averaging
// scheme, that MixingLengthLidDrivenCavitySolver3D already derived and
// validated (see its header for that derivation; it is reused verbatim
// rather than reinvented).
//
// l_s = min(Cs * Delta, kappa * d_wall).
//
// **On the near-wall term, and why it is not van Driest**: plain
// Smagorinsky over-predicts nu_sgs approaching a wall, where the subgrid
// motions it models are suppressed. The textbook remedy is van Driest
// damping, 1 - exp(-y+/A+), which needs y+ -- and therefore the local
// friction velocity u_tau. This project has repeatedly declined to
// estimate u_tau in general geometry (see KEpsilonLidDrivenCavitySolver2D's
// header for the same decision and its reasoning: no exact relation exists
// for a recirculating flow, and estimating it from an unrefined near-wall
// gradient would be exactly the kind of unvalidatable approximation this
// project avoids). Instead the mixing-length cap kappa*d_wall is reused:
// it needs only the exactly-known distance to the nearest of the six
// walls, it vanishes at the wall as required, and it is the same
// self-justifiable geometric bound already validated elsewhere here. This
// is a deliberate, documented substitution -- it is NOT van Driest and is
// not claimed to reproduce it.
//
// **Recalled constant, flagged honestly**: Cs = 0.17 is an empirical
// literature value (commonly quoted between 0.1 and 0.2, and genuinely
// flow-dependent -- one of plain Smagorinsky's well-known weaknesses, and
// the motivation for the dynamic Germano procedure that computes Cs from
// the resolved field instead). It is recalled here, not derived, and is
// exposed as a constructor parameter so it can be varied rather than
// buried. kappa = 0.41 is the same von Karman constant already used
// throughout this project's other closures.
//
// Same ghost-mirror walls, tapered lid, matrix-free CG pressure
// projection, and "no literature benchmark" validation philosophy as every
// other solver here. nu_sgs is exactly zero at all six walls by
// construction, and -- because it is a purely algebraic function of the
// instantaneous strain rate, with no transport equation of its own --
// exactly zero everywhere whenever the velocity field is exactly zero
// (the same strong rest-state guarantee MixingLengthLidDrivenCavitySolver3D
// has, and that the two-equation closures cannot offer).
//
// The staggered-grid plumbing lives in StaggeredCavityBase3D; only the
// subgrid model itself is here.
class SmagorinskyLesLidDrivenCavitySolver3D : public StaggeredCavityBase3D {
public:
    // smagorinskyConstant defaults to the commonly quoted 0.17; see the
    // class comment on why it is a parameter rather than a hardcoded
    // constant.
    SmagorinskyLesLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                           double lengthY, double lengthZ, double viscosity,
                                           double lidVelocity, double smagorinskyConstant = 0.17);

    void step(double dt);

    // The subgrid-scale eddy viscosity nu_sgs at a cell center.
    double subgridViscosity(std::size_t i, std::size_t j, std::size_t k) const;
    // The filter width Delta = cbrt(dx*dy*dz). Constant for a uniform
    // mesh; exposed because the mesh-refinement limit that defines LES is
    // stated in terms of it.
    double filterWidth() const;

private:
    void updateSubgridViscosity();

    double smagorinskyConstant_;
    std::vector<double> nut_;
};

} // namespace aether::solver
