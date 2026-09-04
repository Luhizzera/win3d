#pragma once

#include "aether/solver/StaggeredCavityBase3D.hpp"

#include <cstddef>
#include <vector>

namespace aether::solver {

// Module 6, DES branch: Detached-Eddy Simulation, the classic hybrid
// RANS/LES formulation (Strelets 2001's SST-DES variant of Spalart's
// original 1997 DES), now unblocked because both halves already exist on
// the same 3D staggered cavity: KOmegaSSTLidDrivenCavitySolver3D (RANS,
// near-wall) and SmagorinskyLesLidDrivenCavitySolver3D (LES, away from
// walls).
//
// **The idea in one line**: take the k-omega SST closure verbatim (same
// transport equations, same blending, same near-wall omega condition, same
// Bradshaw-limited eddy-viscosity formula), and replace the RANS
// turbulence length scale L_RANS = sqrt(k)/(beta_star*omega) implicit in
// the k-equation's destruction term with L_DES = min(L_RANS, C_DES*Delta),
// Delta = cbrt(dx*dy*dz) being the exact same LES filter width
// SmagorinskyLesLidDrivenCavitySolver3D already uses.
//
// Concretely: destruction_k = beta_star*k*omega = k^1.5/L_RANS ordinarily.
// Substituting L_DES for L_RANS (algebraically) means multiplying the
// destruction term by F_DES = max(L_RANS/(C_DES*Delta), 1.0):
//   destruction_k = beta_star*k*omega*F_DES.
// F_DES = 1 (plain SST, unmodified) wherever the modeled RANS length scale
// is already smaller than C_DES*Delta -- i.e. wherever the mesh is coarse
// relative to the turbulence it is modeling, which is everywhere on a mesh
// too coarse to resolve any eddies directly. F_DES > 1 wherever the mesh is
// fine enough that C_DES*Delta undercuts L_RANS: destruction grows, k (and
// therefore nu_t, through the same k/omega eddy-viscosity formula k-omega
// SST already uses) drops toward LES-like values. This is exactly the same
// nu_t-shrinks-under-refinement mechanism that defines
// SmagorinskyLesLidDrivenCavitySolver3D, now gated by whichever of L_RANS
// or C_DES*Delta is smaller, instead of applying unconditionally
// everywhere the way plain LES does.
//
// **C_DES = 0.61 is a recalled literature constant, flagged honestly**
// (Strelets' SST-DES value; the original Spalart-Allmaras DES97 paper used
// 0.65 for its own, differently-defined length scale -- not directly
// comparable). Exposed as a constructor parameter, same practice as
// SmagorinskyLesLidDrivenCavitySolver3D's Cs.
//
// **Known, textbook limitation, not fixed here**: this is classic DES97,
// not Delayed DES (DDES). On a mesh where near-wall cells are fine enough
// in the wall-parallel directions while the boundary layer itself is not
// yet resolved, C_DES*Delta can undercut L_RANS *inside* the boundary
// layer, triggering premature LES-mode switching there -- "modeled-stress
// depletion" / grid-induced separation, the well-documented motivation for
// DDES's shielding function. Not built here: this project's 3D grids so
// far are uniform Cartesian (no wall-normal stretching), so the failure
// mode is present in principle but not exercised by the validation used
// (an interior-dominated mesh-refinement comparison against plain SST, not
// a wall-resolved boundary-layer case).
// The staggered-grid plumbing lives in StaggeredCavityBase3D;
// only the closure itself is here.
class DesSstLidDrivenCavitySolver3D : public StaggeredCavityBase3D {
public:
    // cDes defaults to Strelets' commonly quoted 0.61; see the class
    // comment on why it is a parameter rather than a hardcoded constant.
    // useGpu: opt-in (Fase 4 of ROADMAP.md) -- see StaggeredCavityBase3D's
    // own constructor comment for the fallback semantics when GPU is
    // requested but unavailable. The base's momentum predictor already
    // forwards whatever setEddyViscosityField() registered to the GPU
    // path, so this closure's nu_t is honored identically to the CPU one.
    DesSstLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX, double lengthY,
                                   double lengthZ, double viscosity, double lidVelocity, double cDes = 0.61,
                                   bool useGpu = false);

    void step(double dt);

    double k(std::size_t i, std::size_t j, std::size_t k) const;
    double omega(std::size_t i, std::size_t j, std::size_t k) const;
    double eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const;
    // F_DES = max(L_RANS/(C_DES*Delta), 1.0) at a cell center. 1.0 means
    // plain-SST (RANS) behavior there; >1.0 means the k-equation
    // destruction term (and therefore k and nu_t) is being pulled toward
    // LES-like values. Reflects the most recently computed step (0.0
    // everywhere before the first step()).
    double desFactor(std::size_t i, std::size_t j, std::size_t k) const;
    // The LES filter width Delta = cbrt(dx*dy*dz), constant on this uniform
    // mesh; exposed for the same reason
    // SmagorinskyLesLidDrivenCavitySolver3D exposes it.
    double filterWidth() const;

private:
    double kAt(long long i, long long j, long long k) const;
    double omegaGhostAt(std::size_t i, std::size_t j, std::size_t k, int di, int dj, int dk) const;
    double sigmaKNutAt(long long i, long long j, long long k) const;
    double sigmaOmegaNutAt(long long i, long long j, long long k) const;

    void updateBlendingAndCoefficients();

    double cDes_;
    std::vector<double> k_;
    std::vector<double> omega_;
    std::vector<double> nut_;
    std::vector<double> production_;
    std::vector<double> crossDiffusion_;
    std::vector<double> sigmaK_;
    std::vector<double> sigmaOmega_;
    std::vector<double> beta_;
    std::vector<double> gamma_;
    std::vector<double> fDes_;
};

} // namespace aether::solver
