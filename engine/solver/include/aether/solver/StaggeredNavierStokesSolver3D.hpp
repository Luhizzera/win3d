#pragma once

#include <cstddef>
#include <vector>

namespace aether::solver {

// Solves the 3D incompressible Navier-Stokes equations (constant density,
// kinematic viscosity nu) on a triply-periodic **staggered** (Marker-and-
// Cell / MAC) grid via Chorin's projection method -- the same predictor +
// pressure-projection structure as TaylorGreenVortexSolver2D, but now in
// 3D and on a staggered rather than collocated arrangement.
//
// Staggering: pressure sits at cell centers; u, v, w each sit on their own
// grid of face centers (u on x-faces, v on y-faces, w on z-faces), not at
// cell centers. Concretely, for cell (i,j,k) spanning
// [i*dx,(i+1)*dx] x [j*dy,(j+1)*dy] x [k*dz,(k+1)*dz]:
//   u(i,j,k) lives at ( i*dx,       (j+0.5)*dy, (k+0.5)*dz )  -- the cell's
//                                                                 x-minus face
//   v(i,j,k) lives at ( (i+0.5)*dx,  j*dy,      (k+0.5)*dz )  -- y-minus face
//   w(i,j,k) lives at ( (i+0.5)*dx, (j+0.5)*dy,  k*dz      )  -- z-minus face
//   p(i,j,k) lives at ( (i+0.5)*dx, (j+0.5)*dy, (k+0.5)*dz )  -- cell center
// All four fields still use nx*ny*nz entries in a periodic domain (face nx
// coincides with face 0, same as cell nx would).
//
// This is the direct fix for the checkerboard/null-space gap flagged since
// TaylorGreenVortexSolver2D and LidDrivenCavitySolver2D were built: on a
// collocated grid, the compact pressure-gradient and divergence operators
// are each centered *differences* (2h-wide), so their composition does not
// exactly invert the compact Laplacian solved for pressure, leaving a
// residual divergence (~4e-3 measured there). Staggering removes this by
// construction: the pressure gradient at a velocity face uses the two
// *directly adjacent* cell pressures (h-wide, not 2h), and divergence at a
// cell uses the *directly adjacent* velocity faces the same way -- these
// two operators are exact discrete adjoints, so their composition is
// exactly the compact 7-point Laplacian with no null space beyond the
// single constant mode (removed the same way as before, by pinning one
// cell). The convection term needs face-to-face interpolation for the
// cross velocity components (e.g. v and w must be interpolated to a u-face
// location to form the uv/uw flux terms) -- the standard MAC discretization,
// more bookkeeping than the collocated scheme but not conceptually new.
//
// Triply periodic only (no solid walls yet) so it can be validated against
// an exact solution: a 3D Beltrami flow (curl(u) = u), which makes the
// nonlinear advection term collapse to a pure gradient (a standard vector
// identity, (u.grad)u = grad(|u|^2/2) - u x curl(u), and u x u = 0 when
// curl(u) = u) -- so, like TaylorGreenVortexSolver2D's 2D vortex, the whole
// nonlinear 3D velocity field decays as a single exponential in time, an
// exact closed-form solution of the full equations that genuinely exercises
// 3D vortex stretching (unlike embedding a z-invariant 2D flow in a 3D
// grid, which would never touch the z-direction terms).
//
// Solid-wall boundary conditions (a 3D lid-driven cavity) are the natural
// next step once this periodic staggered core is validated, mirroring how
// LidDrivenCavitySolver2D followed TaylorGreenVortexSolver2D.
class StaggeredNavierStokesSolver3D {
public:
    StaggeredNavierStokesSolver3D(std::size_t nx, std::size_t ny, std::size_t nz, double lengthX,
                                   double lengthY, double lengthZ, double viscosity);

    // Sets the velocity components at cell (i,j,k) (used for the initial
    // condition). Per the staggering above, u/v/w are each set at *their
    // own* face location for this (i,j,k), not a shared cell center --
    // callers evaluating an analytic field must sample each component at
    // its own position.
    void setVelocity(std::size_t i, std::size_t j, std::size_t k, double u, double v, double w);

    // A conservative explicit-stability estimate combining the diffusive
    // von Neumann limit with a convective CFL limit based on the given
    // characteristic velocity scale. Callers should use a margin below
    // this (e.g. half), not the exact marginal value.
    double stableTimeStep(double velocityScale) const;

    // Advances one step (predictor + pressure projection) and advances
    // time() by dt.
    void step(double dt);

    double u(std::size_t i, std::size_t j, std::size_t k) const;
    double v(std::size_t i, std::size_t j, std::size_t k) const;
    double w(std::size_t i, std::size_t j, std::size_t k) const;
    double pressure(std::size_t i, std::size_t j, std::size_t k) const;
    double time() const { return time_; }

    // Maximum absolute discrete divergence over all cells for the
    // *current* velocity field. On this staggered grid this should be
    // near machine precision after step() (not just "small and bounded"
    // as on the collocated solvers), since the divergence and pressure-
    // gradient operators are exact discrete adjoints here.
    double maxDivergence() const;

private:
    std::size_t index(std::size_t i, std::size_t j, std::size_t k) const {
        return i + j * nx_ + k * nx_ * ny_;
    }
    std::size_t wrap(long long i, std::size_t n) const;

    std::vector<double> applyPeriodicLaplacian(const std::vector<double>& x) const;
    static double dot(const std::vector<double>& a, const std::vector<double>& b);
    void projectToDivergenceFree(std::vector<double>& uStar, std::vector<double>& vStar,
                                  std::vector<double>& wStar, double dt);

    std::size_t nx_;
    std::size_t ny_;
    std::size_t nz_;
    double dx_;
    double dy_;
    double dz_;
    double viscosity_;
    std::vector<double> u_;
    std::vector<double> v_;
    std::vector<double> w_;
    std::vector<double> p_;
    double time_ = 0.0;
};

} // namespace aether::solver
