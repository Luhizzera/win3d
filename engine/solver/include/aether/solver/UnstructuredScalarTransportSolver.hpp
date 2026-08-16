#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"
#include "aether/solver/UnstructuredFvmBase.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Steady convection-diffusion of a passive scalar on a tetrahedral mesh,
// carried by a **prescribed** velocity field:
//
//   div(U phi) = Gamma laplacian(phi) + S
//
// **Why this class exists, and it is not to simulate anything.** It is the
// ruler for the convection scheme. DIVIDA_TECNICA.md 3.1 records that this
// project's unstructured convection is first-order upwind, chosen for
// stability and known to be the accuracy ceiling of any quantitative result
// -- but "known" was inference, exactly as second order was inference for the
// Laplacian before item 3.2 measured it. Measuring a convection scheme needs
// a case whose exact answer is known, and the Navier-Stokes solver is the
// wrong place to get one: its velocity is what is being solved for, so a
// manufactured solution there measures the momentum equation, the projection
// and the convection scheme all at once.
//
// Here the velocity is an input. Choose any smooth phi, derive the source it
// satisfies, impose phi on the boundary, and the only thing left is the
// discretization -- with the convective and diffusive parts separable by
// turning Gamma down until convection dominates.
//
// **The time march is a means, not the point.** The steady equation with
// convection is not symmetric, so the project's Conjugate Gradient does not
// apply to it directly. Rather than add a non-symmetric solver, this marches
// in time to steady state with convection explicit and diffusion implicit --
// which is exactly how UnstructuredCavitySolver3D treats its momentum
// equation, so the scheme is verified in the form it is actually used, not
// in a form that happens to be easier to solve.
class UnstructuredScalarTransportSolver : public UnstructuredFvmBase {
public:
    // `velocity` is sampled at face centroids to build the face fluxes once.
    // A divergence-free field is expected: the discrete sum of face fluxes
    // around a cell is then zero to O(h^2) rather than exactly, which is a
    // second-order consistency error and so cannot masquerade as a
    // first-order scheme -- the same reasoning the manufactured Laplacian
    // case applies to its midpoint quadrature.
    UnstructuredScalarTransportSolver(const mesh::TetrahedralMesh& mesh, double diffusivity,
                                       std::function<core::Vector3(const core::Vector3&)> velocity,
                                       ConvectionScheme scheme = ConvectionScheme::LimitedLinearUpwind);

    // Boundary faces start insulated; selected ones are held at `value`,
    // evaluated at the face centroid. A convection problem needs the value on
    // *inflow* faces to be well posed; prescribing it on outflow faces too is
    // what a manufactured solution wants, since the exact answer is known
    // there as well.
    void setDirichletBoundary(const std::function<bool(const core::Vector3&)>& selector,
                               const std::function<double(const core::Vector3&)>& value);

    void setSourceTerm(const std::function<double(const core::Vector3&)>& source);

    // The strictest of the per-cell convective limits -- the step a global
    // march would have to take. Reported as a diagnostic; the solve itself
    // uses one step per cell (see localTimeSteps in the implementation).
    double stableTimeStep() const;

    // Marches until the largest per-cell change over a step falls below
    // `tolerance`. Returns the number of steps taken; `lastChange()` says
    // whether it settled or hit the cap, which mean different things.
    std::size_t solveSteady(double tolerance = 1e-12, std::size_t maxSteps = 500000);

    double value(std::size_t cell) const { return phi_.at(cell); }
    double lastChange() const { return lastChange_; }

    // Largest cell Peclet number |U| h / Gamma over the interior faces: how
    // far this case is into the convection-dominated regime, which is the
    // only regime where the convection scheme is what the error measures.
    double maxCellPeclet() const;

private:
    std::vector<double> convectiveFlux(const std::vector<double>& phi) const;
    std::vector<double> localTimeSteps() const;
    std::vector<double> applyTransportOperator(const std::vector<double>& x,
                                                const std::vector<double>& dt) const;

    double diffusivity_;
    std::function<core::Vector3(const core::Vector3&)> velocity_;
    ConvectionScheme scheme_;

    std::vector<bool> boundaryIsDirichlet_;
    std::vector<double> boundaryValue_;
    // U . A through each face, built once: the velocity is prescribed, so
    // these never change.
    std::vector<double> interiorMassFlux_;
    std::vector<double> boundaryMassFlux_;
    std::vector<double> integratedSource_;
    std::vector<double> phi_;
    double lastChange_ = 0.0;
};

} // namespace aether::solver
