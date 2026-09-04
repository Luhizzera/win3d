#pragma once

#include "aether/core/Vector3.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"
#include "aether/solver/UnstructuredFvmBase.hpp"

#include <cstddef>
#include <functional>
#include <vector>

namespace aether::solver {

// Steady diffusion (Laplace's equation) by the finite-volume method on an
// **unstructured tetrahedral mesh** -- the first solver in this project to
// consume mesh generation rather than a uniform Cartesian grid.
//
// **Why this class is the point of ROADMAP Fase 2.** Every other solver
// here takes StructuredGrid3D, so until now the engine could import an STL,
// tetrahedralize it, and then solve in a box. This closes that gap for the
// simplest physics the project has, deliberately: diffusion already has an
// independent analytical answer (the 2D Fourier series that
// SteadyDiffusionSolver is validated against), so solving the *same*
// problem on tetrahedra isolates what changed -- the geometric
// discretization -- from what did not -- the physics.
//
// **What is left in this class is only the physics.** The geometry (the
// over-relaxed decomposition, the non-orthogonal deferred correction, the
// least-squares gradients with their Green-Gauss fallback, the Laplacian and
// the Conjugate Gradient that solves it) is UnstructuredFvmBase, shared with
// the Navier-Stokes solver. What remains here is the boundary condition --
// which faces are held at a fixed value and which are insulated -- and the
// right-hand side that follows from it. That split is the whole reason the
// base exists; see its header.
//
// **The non-orthogonality this rests on.** On a Cartesian grid the face
// normal is parallel to the line joining the two cell centres, so the
// diffusive flux through a face is just (phi_N - phi_P)/|d| times the area.
// On tetrahedra it never is, and that misalignment is a property of mesh
// *shape*, not mesh *size* -- refining a tet mesh does not make it more
// orthogonal (maxNonOrthogonality() stays ~1.5 on these meshes). Fase 2.2
// measured what dropping the resulting correction costs: rms error 2.879 ->
// 2.433 -> 2.365 for n = 4/6/8, an observed order of 0.42 then 0.10 -- a
// plateau, not convergence. Running that experiment before committing to the
// simpler scheme is the Fase 1 lesson applied.
//
// **With the correction in place this class is second order, and that is now
// measured rather than inferred.** Fase 2.2 could only report 0.91 -> 0.98 on
// the plate problem, whose hot and cold edges meet at two corners where the
// exact solution has no bounded gradient -- no scheme converges at its formal
// order in a norm that includes those cells, so the plate could not tell a
// second-order discretization from a first-order one. A manufactured solution
// smooth over the whole closed cube gives order 2.013, 2.024, 2.026 for
// n = 4/6/8/10 (and 1.947 at n = 12), on meshes whose non-orthogonality stays
// between 1.45 and 1.75. See testManufacturedSolutionReachesSecondOrder() and
// setSourceTerm(), which exists for it.
class UnstructuredDiffusionSolver : public UnstructuredFvmBase {
public:
    explicit UnstructuredDiffusionSolver(const mesh::TetrahedralMesh& mesh);

    // Every boundary face starts insulated (zero-gradient / Neumann), which
    // contributes no flux and therefore needs no coefficient at all. Faces
    // whose *centroid* satisfies `selector` are instead fixed to `value`.
    // Selecting by centroid rather than by index is what lets a caller say
    // "the face on the x = 0 plane" without knowing anything about how the
    // mesh generator happened to number things. A later call overrides an
    // earlier one for the same face.
    void setDirichletBoundary(const std::function<bool(const core::Vector3&)>& selector, double value);

    // The same, with the value varying over the boundary: `value` is
    // evaluated at each selected face's centroid. Needed by any problem whose
    // boundary data is not piecewise constant -- a manufactured solution
    // above all, where the boundary condition *is* the exact solution
    // sampled on the surface. The constant overload is this one with a
    // constant function, and is kept because "this whole side is held at
    // 100 degrees" should not have to be written as a lambda.
    void setDirichletBoundary(const std::function<bool(const core::Vector3&)>& selector,
                               const std::function<double(const core::Vector3&)>& value);

    // Volumetric source: solves the Poisson equation
    //
    //   laplacian(phi) + S(x) = 0
    //
    // instead of Laplace's. `source` is evaluated at each cell centroid and
    // multiplied by the cell volume, which is the midpoint rule for the
    // integral the finite-volume method actually needs. That rule is exact
    // for a linear source and carries an O(h^2) relative error otherwise --
    // deliberately noted, because it sits at exactly the order this class's
    // accuracy is claimed at, so a coarser quadrature here would cap the
    // measured convergence rate and look like a defect of the scheme.
    //
    // **Added for the method of manufactured solutions** (ROADMAP Fase 2.2
    // left second order as an inference, DIVIDA_TECNICA.md 3.2): pick any
    // smooth phi, set S = -laplacian(phi), impose phi on the boundary, and
    // the discretization error is the only thing left to measure. Without a
    // source term the only exactly-known solutions are harmonic ones, which
    // is a narrower family and leaves the volume integration untested.
    //
    // Pass {} to clear it and return to Laplace's equation.
    void setSourceTerm(std::function<double(const core::Vector3&)> source);

    // Per-cell thermal conductivity, evaluated once at each cell centroid --
    // same convention as setSourceTerm. Default: uniform 1.0, which
    // reproduces solveConjugateGradient() bit-for-bit and is what makes this
    // opt-in. A "solid" region and a "fluid" region are just two different
    // values of k selected by position, the same convention wallVelocity/
    // isOutlet/setDirichletBoundary already use.
    //
    // **This is ROADMAP Fase 6.1 (conjugate heat transfer) in its literal
    // scope: solid and fluid regions with different conductivities, coupled
    // at their shared interface.** The interior face conductivity is the
    // distance-weighted harmonic mean of the two neighbouring cells' values
    // -- the two-point-equivalent coefficient that continuity of flux across
    // the interface implies, derived in updateConductivity()'s own comment.
    // A boundary face uses its owning cell's conductivity, since a
    // prescribed value sits exactly at the face with nothing on the other
    // side. Deliberately conduction only, no velocity anywhere: this isolates
    // the one thing that is new (the interface coupling) from advection,
    // which UnstructuredCavitySolver3D's EnergyModel already covers.
    // Integrating this into a real flow around a solid obstacle is a
    // separate, natural extension, not required by this scope.
    //
    // Pass {} to clear it and return to uniform conductivity 1.0.
    void setConductivity(std::function<double(const core::Vector3&)> conductivity);

    // Outer deferred-correction sweeps, each an inner matrix-free CG solve
    // (same pattern as SteadyDiffusionSolver). Returns the number of outer
    // sweeps actually run; sweeps stop early once the solution stops moving
    // by more than `tolerance` between them.
    //
    // Requires at least one Dirichlet face: with pure Neumann boundaries the
    // operator is singular (the solution is defined only up to a constant),
    // and this class does not pin a reference cell the way the pressure
    // solvers do -- the diffusion problems it is meant for always have a
    // fixed value somewhere.
    std::size_t solveConjugateGradient(std::size_t maxIterations = 20000, double tolerance = 1e-10,
                                        std::size_t maxOuterSweeps = 50);

    // Gradient at each cell from the current solution. Exposed because it is
    // the quantity the non-orthogonal correction is built from, so a caller
    // checking this class's accuracy can inspect it.
    std::vector<core::Vector3> cellGradients() const { return computeCellGradients(phi_); }

    // Largest per-cell change across the final outer sweep. The number that
    // says whether the deferred-correction loop actually converged or merely
    // ran out of sweeps -- without it, a solve that stopped at the iteration
    // cap is indistinguishable from one that settled, and the two mean very
    // different things about the result.
    double lastOuterChange() const { return lastOuterChange_; }

    double value(std::size_t cell) const { return phi_.at(cell); }

private:
    void rebuildCoefficients();

    // Rebuilds interiorFaceWeight_/boundaryFaceWeight_ from cellConductivity_
    // -- see setConductivity()'s declaration for the physics. Called from the
    // constructor (uniform 1.0) and from setConductivity().
    void updateConductivity();

    // Per-boundary-face prescribed value, or "no value" when insulated.
    // Indexed by *mesh* face, since that is what a selector resolves to and
    // what the base hands back in BoundaryFace::meshFace.
    std::vector<bool> boundaryIsDirichlet_;
    std::vector<double> boundaryValue_;
    bool hasDirichletFace_ = false;
    // Already integrated: S(centroid) * V, one per cell, so the solve does
    // not re-evaluate the caller's function on every outer sweep.
    std::vector<double> integratedSource_;
    std::vector<double> phi_;
    double lastOuterChange_ = 0.0;

    // Only true after setConductivity() is called with a non-empty function:
    // while false, solveConjugateGradient() runs the original unweighted
    // path (applyLaplacian/solveDeferredCorrection) exactly as before this
    // capability existed, so no existing caller's result changes at all.
    bool hasConductivity_ = false;
    std::vector<double> cellConductivity_;      // one per cell; uniform 1.0 by default
    std::vector<double> interiorFaceWeight_;    // parallel to interiorFaces_
    std::vector<double> boundaryFaceWeight_;    // parallel to boundaryFaces_
};

} // namespace aether::solver
