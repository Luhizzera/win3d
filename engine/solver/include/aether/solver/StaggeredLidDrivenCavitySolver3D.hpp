#pragma once

#include "aether/solver/StaggeredCavityBase3D.hpp"

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
// The staggered-grid plumbing (indexing, ghost mirrors, momentum
// predictor, pressure projection, stableTimeStep, maxDivergence) lives in
// StaggeredCavityBase3D. Nothing else is left here: with no turbulence
// closure registered, nutAt() is identically zero and the base's
// face-weighted diffusion reduces to the plain molecular Laplacian this
// class used to spell out itself.
class StaggeredLidDrivenCavitySolver3D : public StaggeredCavityBase3D {
public:
    // `convection` defaults to Central, reproducing every number this class
    // has ever produced bit-for-bit -- see DIVIDA_TECNICA.md 4.4 for why,
    // and StaggeredCavityBase3D::ConvectionScheme for what the other two
    // options mean on this staggered grid specifically. Exposed here (and
    // not yet on the five turbulence closures) because this is the class
    // the Re=400 measurement that would justify changing the default runs
    // against -- the same one 2D's own port was measured on.
    // useGpu: opt-in (Fase 4 of ROADMAP.md) -- see
    // StaggeredCavityBase3D's own constructor comment for the fallback
    // semantics when GPU is requested but unavailable. Exposed on this
    // class specifically for the same reason `convection` already is: this
    // is the one the GPU speedup was measured against.
    StaggeredLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX, double lengthY, double lengthZ,
           double viscosity, double lidVelocity,
           ConvectionScheme convection = ConvectionScheme::Central, bool useGpu = false);

    void step(double dt);

    // Safe to expose as-is (unlike the five turbulent closures): this
    // class carries no state beyond what StaggeredCavityBase3D::loadState
    // already covers (u, v, w, p, time) -- no eddy-viscosity field or other
    // closure-specific state to also restore.
    using StaggeredCavityBase3D::loadState;
};

} // namespace aether::solver
