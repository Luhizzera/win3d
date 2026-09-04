#pragma once

#include "aether/solver/StaggeredCavityBase3D.hpp"

#include <cstddef>
#include <vector>

namespace aether::solver {

// k-omega SST extended from the real 2D cavity (KOmegaSSTLidDrivenCavitySolver2D)
// to the real 3D solid-wall cavity, mirroring the same step
// KEpsilonLidDrivenCavitySolver3D just took for the simpler two-equation
// closure. Same self-derived omega wall condition as the 2D class
// (omega_wall = 2*nu/(beta**y_half^2), from substituting the general
// near-wall k~y^2 relation into the exact epsilon=beta**k*omega identity
// and taking y->0 -- see KOmegaSSTLidDrivenCavitySolver2D's own comment for
// the full derivation and the acknowledged departure from Wilcox's more
// careful matched-asymptotic formula), same blended F1/F2 SST formulation,
// same Patankar-linearized destruction terms, same warm-start practice
// (now from MixingLengthLidDrivenCavitySolver3D).
//
// **The one genuinely new physics 3D introduces: vorticity is a vector,
// not a scalar.** The 2D class's Bradshaw shear-stress limiter uses
// |dv/dx - du/dy|, the z-component of curl(u) -- the *only* nonzero
// component in a genuinely 2D flow. In 3D, curl(u) = (dw/dy-dv/dz,
// du/dz-dw/dx, dv/dx-du/dy) has all three components in general, so the
// limiter uses the full vorticity magnitude sqrt(wx^2+wy^2+wz^2). Computed
// from the same 6 cross-derivatives (dudy, dvdx, dudz, dwdx, dvdz, dwdy)
// already needed for the strain-rate tensor -- reusing
// MixingLengthLidDrivenCavitySolver3D's edge-averaging idiom once again,
// just also combining the same six numbers antisymmetrically (vorticity)
// instead of only symmetrically (strain).
//
// k, omega, nu_t, production, cross-diffusion and the blended coefficients
// all live at cell centers (same as KEpsilonLidDrivenCavitySolver3D), so
// their transport uses the same plain colocated 6-neighbor stencil; only
// the momentum equations need the staggered-grid own-axis/transverse
// treatment.
// The staggered-grid plumbing (indexing, ghost mirrors, momentum
// predictor, pressure projection, stableTimeStep, maxDivergence) lives in
// StaggeredCavityBase3D; only the closure itself is here.
class KOmegaSSTLidDrivenCavitySolver3D : public StaggeredCavityBase3D {
public:
    // useGpu: opt-in (Fase 4 of ROADMAP.md) -- see StaggeredCavityBase3D's
    // own constructor comment for the fallback semantics when GPU is
    // requested but unavailable. The base's momentum predictor already
    // forwards whatever setEddyViscosityField() registered to the GPU
    // path, so this closure's nu_t is honored identically to the CPU one.
    KOmegaSSTLidDrivenCavitySolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                      double lengthY, double lengthZ, double viscosity, double lidVelocity,
                                      bool useGpu = false);

    void step(double dt);

    double k(std::size_t i, std::size_t j, std::size_t k) const;
    double omega(std::size_t i, std::size_t j, std::size_t k) const;
    double eddyViscosity(std::size_t i, std::size_t j, std::size_t k) const;

private:
    // Cell-centered k; homogeneous Dirichlet ghost mirror (k=0 at any wall).
    double kAt(long long i, long long j, long long k) const;
    // omega's ghost mirrors the computed wall value 2*nu/(beta**y_half^2).
    double omegaGhostAt(std::size_t i, std::size_t j, std::size_t k, int di, int dj, int dk) const;
    // sigma_k(or omega) * nu_t at a neighbor cell; 0.0 outside the domain.
    double sigmaKNutAt(long long i, long long j, long long k) const;
    double sigmaOmegaNutAt(long long i, long long j, long long k) const;

    void updateBlendingAndCoefficients();

    std::vector<double> k_;
    std::vector<double> omega_;
    std::vector<double> nut_;
    std::vector<double> production_;
    std::vector<double> crossDiffusion_;
    std::vector<double> sigmaK_;
    std::vector<double> sigmaOmega_;
    std::vector<double> beta_;
    std::vector<double> gamma_;
};

} // namespace aether::solver
