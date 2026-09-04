#pragma once

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Interface tracking by the level-set method -- ROADMAP Fase 6.3, and
// deliberately only the first of the two things that phase names
// ("rastreamento de interface, tensão superficial"). No surface tension,
// no coupling to a solved two-phase Navier-Stokes: the velocity field this
// class advects an interface through is prescribed by the caller, not
// computed. Full multiphase flow is comparable in scope to this project's
// entire engine, not one step -- this class is the smallest genuinely new
// piece: an interface that keeps a correct shape and position under a
// known flow, which plain scalar advection elsewhere in this project has
// never actually been asked to demonstrate (a translating/rotating shape
// returning to itself, not just a passive field settling into some
// profile).
//
// **Level set, not volume-of-fluid.** VOF's conservative form
// d(C)/dt + div(u*C) = 0 needs a geometric interface reconstruction (PLIC)
// or a compressive scheme to keep from smearing. Level set's
// non-conservative form dphi/dt + u.grad(phi) = 0 needs neither -- just an
// advection equation, at the cost of needing periodic reinitialization to
// keep phi close to a true signed distance (|grad(phi)| = 1) over long
// runs. Chosen here because it needs no new geometric algorithm, and
// because reinitialization only affects the field *away* from the
// interface, not the interface's own position under pure advection -- so
// this first slice, validated over exactly one rotation period, does not
// need it yet. Historically the simplest member of this family (Osher &
// Sethian 1988), the same reasoning that already chose Rusanov before
// Roe/HLLC and mixing length before k-epsilon elsewhere in this project.
//
// **Numerical scheme: first-order upwind finite difference**, the
// classical Hamilton-Jacobi discretization -- there is no face flux to
// integrate (the equation is non-conservative), so this does not reuse
// ConvectionLimiter.hpp's face-value blend, which is shaped for exactly
// that (phi_f = phi_C + psi*(phi_central-phi_C)) and does not apply to a
// bare derivative. A one-sided difference chosen by the local sign of u
// (respectively v) at each cell, forward Euler in time -- first order in
// both, deliberately: the simplest pair that is still correct, matching
// this project's own precedent (Rusanov+forward-Euler for Fase 6.2).
// Second-order (ENO/WENO) is a natural extension, not this class's scope.
class LevelSetAdvectionSolver2D {
public:
    LevelSetAdvectionSolver2D(std::size_t nx, std::size_t ny, double lengthX, double lengthY);

    // Evaluated once per cell centre. Positive inside the tracked region,
    // negative outside -- a genuine signed distance is the natural choice
    // (e.g. radius minus distance to a circle's centre) but not enforced;
    // this class advects whatever scalar it is given.
    void initialize(const std::function<double(double x, double y)>& signedDistance);

    // Prescribed velocity, evaluated once per cell centre at call time --
    // same convention as UnstructuredDiffusionSolver::setSourceTerm/
    // setConductivity (evaluated once, never inside step()). Deliberately
    // not solved from anything: this class tracks an interface through a
    // *given* flow. Coupling to a solved two-phase velocity field is
    // future work, not this class's scope.
    void setVelocityField(const std::function<void(double x, double y, double& u, double& v)>& velocity);

    // One forward-Euler step of the upwind Hamilton-Jacobi update.
    void step(double dt);

    // CFL-limited step: cfl * min(dx,dy) / max_i(|u_i|, |v_i|). Measured on
    // the solid-body rotation validation case: stable through cfl=0.70,
    // diverges (|phi| reaching 1e24 within one period) at cfl=0.75 --
    // *not* the naive cfl<=1 a single-direction bound would suggest.
    // Updating both directions in one un-split step needs the combined
    // bound dt*(|u|/dx+|v|/dy)<=1, which for comparable |u|,|v| and
    // dx==dy allows only about half of what this function's single-
    // direction formula implies at cfl=1 -- consistent with the measured
    // break landing near 0.7-0.75 rather than 1.0. The default keeps
    // roughly 1.75x margin below the measured break.
    double stableTimeStep(double cfl = 0.4) const;

    double value(std::size_t i, std::size_t j) const { return phi_.at(index(i, j)); }
    double cellCenterX(std::size_t i) const;
    double cellCenterY(std::size_t j) const;
    double time() const { return time_; }

    // Area of the region where value() > 0, by a plain cell-count*dx*dy
    // sum. Exact (up to the scheme's own numerical diffusion, measured and
    // reported, never assumed zero) under any divergence-free velocity
    // field, by the Reynolds transport theorem -- the conservation half of
    // this class's validation gate.
    double insideArea() const;

private:
    std::size_t index(std::size_t i, std::size_t j) const { return j * nx_ + i; }

    // Ghost-aware lookup: constant extrapolation (mirrors the nearest real
    // cell unchanged) at all four edges. Sufficient because this class's
    // own validation case never lets the tracked interface approach a
    // boundary, so the treatment there is never actually exercised.
    double phiAt(long long i, long long j) const;

    std::size_t nx_;
    std::size_t ny_;
    double lengthX_;
    double lengthY_;
    double dx_;
    double dy_;
    double time_ = 0.0;
    std::vector<double> phi_;
    std::vector<double> u_;
    std::vector<double> v_;
};

} // namespace aether::solver
