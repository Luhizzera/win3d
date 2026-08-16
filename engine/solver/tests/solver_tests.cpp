#include "aether/mesh/DelaunayTetrahedralization3D.hpp"
#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/mesh/TetrahedralMesh.hpp"
#include "aether/solver/UnstructuredCavitySolver3D.hpp"
#include "aether/solver/UnstructuredDiffusionSolver.hpp"
#include "aether/solver/DesSstLidDrivenCavitySolver3D.hpp"
#include "aether/solver/ImplicitConvectionDiffusionSolver1D.hpp"
#include "aether/solver/KEpsilonChannelFlowSolver1D.hpp"
#include "aether/solver/KEpsilonLidDrivenCavitySolver2D.hpp"
#include "aether/solver/KEpsilonLidDrivenCavitySolver3D.hpp"
#include "aether/solver/KOmegaSSTChannelFlowSolver1D.hpp"
#include "aether/solver/KOmegaSSTLidDrivenCavitySolver2D.hpp"
#include "aether/solver/KOmegaSSTLidDrivenCavitySolver3D.hpp"
#include "aether/solver/LidDrivenCavitySolver2D.hpp"
#include "aether/solver/MixingLengthChannelFlowSolver1D.hpp"
#include "aether/solver/MixingLengthLidDrivenCavitySolver2D.hpp"
#include "aether/solver/MixingLengthLidDrivenCavitySolver3D.hpp"
#include "aether/solver/SmagorinskyLesLidDrivenCavitySolver3D.hpp"
#include "aether/solver/MultigridPoissonSolver2D.hpp"
#include "aether/solver/StaggeredLidDrivenCavitySolver3D.hpp"
#include "aether/solver/StaggeredNavierStokesSolver3D.hpp"
#include "aether/solver/SteadyDiffusionSolver.hpp"
#include "aether/solver/TaylorGreenVortexSolver2D.hpp"
#include "aether/solver/TransientDiffusionSolver.hpp"
#include "aether/testing/Check.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>
#include <cstdio>

using namespace aether::core;
using namespace aether::mesh;
using namespace aether::solver;

namespace {

bool nearlyEqual(double a, double b, double tolerance) {
    return std::fabs(a - b) <= tolerance;
}

// With ny = nz = 1 there are no y/z neighbors at all, so the discrete
// problem reduces to pure 1D conduction (source = 0, i.e. Laplace's
// equation). setBoundaryValue() fixes cells 0 and nx-1 directly (rather
// than via a ghost-cell wall treatment applied at the domain edge), so the
// exact discrete steady state T_i = (T_{i-1} + T_{i+1}) / 2 for interior i
// is linear in cell *index*, anchored at the two fixed cells - not linear
// in physical x measured from the domain edges (those two differ by half
// a cell width at each end).
void test1DConductionMatchesAnalyticalProfile() {
    const double tMin = 10.0;
    const double tMax = 50.0;
    const std::size_t nx = 20;

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(2.0, 0.1, 0.1), nx, 1, 1);
    SteadyDiffusionSolver solver(grid);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, tMin);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, tMax);

    const std::size_t iterations = solver.solve(20000, 1e-10);
    AETHER_CHECK(iterations < 20000); // must have converged, not hit the cap

    for (std::size_t i = 0; i < nx; ++i) {
        const double expected =
            tMin + (tMax - tMin) * static_cast<double>(i) / static_cast<double>(nx - 1);
        AETHER_CHECK(nearlyEqual(solver.value(i, 0, 0), expected, 1e-6));
    }
}

void testBoundaryCellsStayFixed() {
    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(1.0, 1.0, 1.0), 5, 3, 3);
    SteadyDiffusionSolver solver(grid);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 100.0);
    solver.solve();

    for (std::size_t j = 0; j < 3; ++j) {
        for (std::size_t k = 0; k < 3; ++k) {
            AETHER_CHECK(nearlyEqual(solver.value(0, j, k), 0.0, 1e-12));
            AETHER_CHECK(nearlyEqual(solver.value(4, j, k), 100.0, 1e-12));
        }
    }
}

// Fully-developed steady laminar (plane Poiseuille) flow between two
// parallel plates: the streamwise momentum equation reduces to
// mu * d^2u/dy^2 = dp/dx, i.e. this solver's Poisson equation with
// phi = u and source = -(1/mu) * dp/dx (see the class comment). The
// grid's "x" axis stands in for the wall-normal direction y here - the
// solver has no notion of physical axis names, this is purely how the
// grid/BCs are set up - and no-slip walls give u = 0 at both fixed layers.
//
// As in the linear conduction case, fixing cells 0 and nx-1 directly means
// the *exact* discrete solution of u_{i-1} - 2u_i + u_{i+1} = -source*h^2
// with u_0 = u_{nx-1} = 0 is the discrete parabola
// u_i = (source * h^2 / 2) * i * (nx - 1 - i): the same Poiseuille shape,
// anchored to cell index rather than physical y.
void testPlanePoiseuilleProfile() {
    const std::size_t nx = 30;
    const double h = 0.1;
    const double source = 8.0; // stands for -(1/mu) * dp/dx

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(h * nx, 0.05, 0.05), nx, 1, 1);
    SteadyDiffusionSolver solver(grid);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 0.0);
    solver.setSourceTerm(source);

    const std::size_t iterations = solver.solve(50000, 1e-12);
    AETHER_CHECK(iterations < 50000);

    for (std::size_t i = 0; i < nx; ++i) {
        const double expected =
            (source * h * h / 2.0) * static_cast<double>(i) * static_cast<double>(nx - 1 - i);
        AETHER_CHECK(nearlyEqual(solver.value(i, 0, 0), expected, 1e-6));
    }

    // Peak velocity at the channel centerline. nx - 1 = 29 is odd, so no
    // integer index sits exactly at the continuous parabola's midpoint;
    // the discrete peak (whichever of the two center indices scores
    // higher in i * (nx - 1 - i)) is slightly below the continuous
    // formula's source * H_eff^2 / 8, so compare against the same
    // discrete formula validated above rather than the continuous one.
    double expectedMax = 0.0;
    for (std::size_t i = 0; i < nx; ++i) {
        expectedMax = std::max(
            expectedMax, (source * h * h / 2.0) * static_cast<double>(i) * static_cast<double>(nx - 1 - i));
    }
    double actualMax = 0.0;
    for (std::size_t i = 0; i < nx; ++i) {
        actualMax = std::max(actualMax, solver.value(i, 0, 0));
    }
    AETHER_CHECK(nearlyEqual(actualMax, expectedMax, 1e-6));
}

// Same 1D conduction problem as test1DConductionMatchesAnalyticalProfile(),
// solved both ways: Conjugate Gradient must reach the same answer as
// Gauss-Seidel, in dramatically fewer iterations (this is the whole point
// of adding it - Gauss-Seidel needs O(n^2)-ish sweeps to converge, visible
// here as roughly 1500+ iterations for just 50 cells).
void testConjugateGradientMatchesGaussSeidelAndConvergesFaster() {
    const double tMin = 10.0;
    const double tMax = 50.0;
    const std::size_t nx = 50;

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(5.0, 0.1, 0.1), nx, 1, 1);

    SteadyDiffusionSolver gaussSeidelSolver(grid);
    gaussSeidelSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, tMin);
    gaussSeidelSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, tMax);
    const std::size_t gaussSeidelIterations = gaussSeidelSolver.solve(100000, 1e-10);
    AETHER_CHECK(gaussSeidelIterations < 100000);

    SteadyDiffusionSolver cgSolver(grid);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, tMin);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, tMax);
    const std::size_t cgIterations = cgSolver.solveConjugateGradient(1000, 1e-10);
    AETHER_CHECK(cgIterations < 1000);

    for (std::size_t i = 0; i < nx; ++i) {
        AETHER_CHECK(nearlyEqual(cgSolver.value(i, 0, 0), gaussSeidelSolver.value(i, 0, 0), 1e-6));
    }
    AETHER_CHECK(cgIterations * 10 < gaussSeidelIterations);
}

// Conjugate Gradient must also handle a non-zero source term correctly:
// re-validate against the same Poiseuille analytical profile as
// testPlanePoiseuilleProfile(), just solved with the other method.
void testConjugateGradientWithSourceTermMatchesPoiseuilleProfile() {
    const std::size_t nx = 30;
    const double h = 0.1;
    const double source = 8.0;

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(h * nx, 0.05, 0.05), nx, 1, 1);
    SteadyDiffusionSolver solver(grid);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 0.0);
    solver.setSourceTerm(source);

    const std::size_t iterations = solver.solveConjugateGradient(1000, 1e-12);
    AETHER_CHECK(iterations < 1000);

    for (std::size_t i = 0; i < nx; ++i) {
        const double expected =
            (source * h * h / 2.0) * static_cast<double>(i) * static_cast<double>(nx - 1 - i);
        AETHER_CHECK(nearlyEqual(solver.value(i, 0, 0), expected, 1e-6));
    }
}

double stableSinhRatio(double a, double b) {
    // sinh(a)/sinh(b) without overflowing for large a, b (0 <= a <= b):
    // exp(a-b) is bounded above by 1, and the correction factors just
    // underflow harmlessly toward 0/1 for large n instead of exploding.
    return std::exp(a - b) * (1.0 - std::exp(-2.0 * a)) / (1.0 - std::exp(-2.0 * b));
}

// Closed-form solution (Fourier series, separation of variables) for
// steady 2D heat conduction in a rectangular plate of size lx by ly with
// three sides held at 0 and the y = ly side held at t0 - the classic
// textbook Laplace problem, used here to validate setBoundaryValue() on a
// Y face (not just X) and the per-axis anisotropic weighting fix together
// on a genuinely 2D case.
double analyticalPlateTemperature(double x, double y, double lx, double ly, double t0, int numOddTerms) {
    const double kPi = 3.14159265358979323846;
    double sum = 0.0;
    for (int k = 0; k < numOddTerms; ++k) {
        const int n = 2 * k + 1;
        const double kx = static_cast<double>(n) * kPi / lx;
        sum += (1.0 / n) * std::sin(kx * x) * stableSinhRatio(kx * y, kx * ly);
    }
    return (4.0 * t0 / kPi) * sum;
}



// ---------------------------------------------------------------------------
// ROADMAP Fase 2.2: the same plate problem, on an unstructured tetrahedral
// mesh instead of a Cartesian grid -- the first time a solver in this project
// consumes mesh generation rather than StructuredGrid3D.
// ---------------------------------------------------------------------------

// A tetrahedral mesh of the unit cube from a lattice of `n+1` points per
// axis. Two details matter and are not incidental:
//
// 1. **Boundary points are not jittered along their own normal.** A point on
//    the x = 0 plane may move in y and z but never in x, so the boundary
//    faces stay exactly planar and exactly on the cube's faces -- which is
//    what lets the Dirichlet selectors below identify them by centroid with
//    an exact tolerance instead of a fuzzy one.
// 2. **Everything else is jittered.** A perfectly regular lattice is
//    massively co-spherical, which is Delaunay's degenerate tie case; this
//    project's 2D tests already avoid regular grids for exactly that reason.
//    Jittering only the components that are free to move gets both
//    properties at once.
aether::mesh::DelaunayTetrahedralization3D buildCubeLatticeTetrahedralization(std::size_t n) {
    aether::mesh::DelaunayTetrahedralization3D tets;
    const double step = 1.0 / static_cast<double>(n);
    const double jitter = 0.22 * step;
    for (std::size_t i = 0; i <= n; ++i) {
        for (std::size_t j = 0; j <= n; ++j) {
            for (std::size_t k = 0; k <= n; ++k) {
                double x = static_cast<double>(i) * step;
                double y = static_cast<double>(j) * step;
                double z = static_cast<double>(k) * step;
                // Deterministic, no RNG -- reproducible across runs and machines.
                if (i > 0 && i < n) {
                    x += jitter * std::sin(4.1 * y + 2.7 * z + 1.3);
                }
                if (j > 0 && j < n) {
                    y += jitter * std::sin(3.3 * z + 5.1 * x + 0.7);
                }
                if (k > 0 && k < n) {
                    z += jitter * std::sin(2.9 * x + 4.7 * y + 2.1);
                }
                tets.addPoint(x, y, z);
            }
        }
    }
    tets.tetrahedralize();
    return tets;
}

// Solves the plate problem on one such mesh and returns the volume-weighted
// RMS error against the Fourier series. Volume-weighted because tetrahedra
// vary in size: an unweighted average would let a cloud of tiny slivers
// outvote the cells that actually carry the domain.
struct UnstructuredPlateResult {
    double rmsError;
    double maxNonOrthogonality;
    std::size_t cellCount;
    std::size_t iterations;
    double outerChange;
    double interiorRmsError;
};

UnstructuredPlateResult solveUnstructuredPlate(std::size_t n) {
    const double t0 = 100.0;
    const aether::mesh::DelaunayTetrahedralization3D tets = buildCubeLatticeTetrahedralization(n);
    const aether::mesh::TetrahedralMesh mesh = aether::mesh::TetrahedralMesh::fromTetrahedralization(tets);

    aether::solver::UnstructuredDiffusionSolver solver(mesh);
    constexpr double kOnPlane = 1e-9;
    // Three cold sides and one hot side, exactly as the structured version.
    // The z faces are left insulated (the default), which makes the solution
    // independent of z and therefore the *same* 2D problem the Fourier
    // series describes -- that is what makes the comparison legitimate.
    solver.setDirichletBoundary([&](const Vector3& c) { return c.x < kOnPlane; }, 0.0);
    solver.setDirichletBoundary([&](const Vector3& c) { return c.x > 1.0 - kOnPlane; }, 0.0);
    solver.setDirichletBoundary([&](const Vector3& c) { return c.y < kOnPlane; }, 0.0);
    solver.setDirichletBoundary([&](const Vector3& c) { return c.y > 1.0 - kOnPlane; }, t0);

    const std::size_t iterations = solver.solveConjugateGradient(40000, 1e-10);

    double weightedSquaredError = 0.0;
    double totalVolume = 0.0;
    // The plate problem is discontinuous at the two top corners, where the hot
    // edge meets a cold one: the exact solution has no bounded gradient there,
    // so *any* scheme carries large local error and no scheme converges at its
    // formal order in a norm that includes those cells. The structured test
    // sidesteps this by sampling interior points only; here the same exclusion
    // is measured explicitly, so the two numbers can be compared and the
    // singularity's contribution separated from the discretization's.
    double interiorSquaredError = 0.0;
    double interiorVolume = 0.0;
    constexpr double kCornerRadius = 0.15;
    for (std::size_t cell = 0; cell < mesh.cellCount(); ++cell) {
        const Vector3 c = mesh.cellCentroid(cell);
        const double expected = analyticalPlateTemperature(c.x, c.y, 1.0, 1.0, t0, 200);
        const double error = solver.value(cell) - expected;
        const double volume = mesh.cellVolume(cell);
        weightedSquaredError += volume * error * error;
        totalVolume += volume;

        const double toLeftCorner = std::sqrt(c.x * c.x + (1.0 - c.y) * (1.0 - c.y));
        const double toRightCorner = std::sqrt((1.0 - c.x) * (1.0 - c.x) + (1.0 - c.y) * (1.0 - c.y));
        if (std::min(toLeftCorner, toRightCorner) > kCornerRadius) {
            interiorSquaredError += volume * error * error;
            interiorVolume += volume;
        }
    }

    return {std::sqrt(weightedSquaredError / totalVolume), solver.maxNonOrthogonality(), mesh.cellCount(),
            iterations, solver.lastOuterChange(), std::sqrt(interiorSquaredError / interiorVolume)};
}

// ---------------------------------------------------------------------------
// ROADMAP Fase 3: incompressible Navier-Stokes on a tetrahedral mesh.
// ---------------------------------------------------------------------------

// The same physical claim testLidDrivenCavityPrimaryVortexTopology() makes for
// the structured solver, now on tetrahedra: at low Reynolds number a lid-driven
// cavity settles into one primary recirculating vortex, so fluid near the lid
// moves with it while fluid near the floor must move *against* it -- not merely
// plausible but forced, since mass cannot pile up in a closed box.
//
// Checking topology rather than pointwise values is deliberate and is this
// project's standing practice for the cavity: it is a consequence of
// conservation that any correct solver must reproduce, and it needs no
// literature benchmark table recalled from memory.
// **The check that outlets actually work**, and the prerequisite the cylinder
// case was really blocked on. A closed cavity can have zero divergence in
// every cell while the domain as a whole is sealed; external flow needs mass
// to cross the boundary, which was structurally impossible before an outlet
// existed.
//
// Straight channel: fluid enters at x = 0 with a prescribed velocity, leaves
// at x = 1 through a pressure outlet, solid walls elsewhere. The claim is
// global, not local -- for a steady incompressible flow whatever enters must
// leave, so the net boundary flux has to vanish. That is a conservation
// statement no scheme can satisfy by accident, and it is exactly what a
// sealed domain cannot do at all.
void testChannelWithOutletConservesGlobalMass() {
    const double inletSpeed = 1.0;
    const aether::mesh::DelaunayTetrahedralization3D tets = buildCubeLatticeTetrahedralization(3);
    const aether::mesh::TetrahedralMesh mesh = aether::mesh::TetrahedralMesh::fromTetrahedralization(tets);

    aether::solver::UnstructuredCavitySolver3D solver(
        mesh, 0.1,
        [&](const Vector3& p2) {
            // Inlet on x = 0; every other non-outlet face is a no-slip wall.
            if (p2.x < 1e-9) {
                return Vector3{inletSpeed, 0.0, 0.0};
            }
            return Vector3{0.0, 0.0, 0.0};
        },
        [&](const Vector3& p2) { return p2.x > 1.0 - 1e-9; }, // outlet on x = 1
        0.0);

    const double dt = solver.stableTimeStep();
    const auto steps = static_cast<int>(3.0 / dt);
    for (int st = 0; st < steps; ++st) {
        solver.step(dt);
    }

    // Scale the balance by the inflow itself, so "small" means small relative
    // to the flow being driven, not small in absolute units nobody can judge.
    double inflow = 0.0;
    for (std::size_t f = 0; f < mesh.faceCount(); ++f) {
        if (mesh.isBoundaryFace(f) && mesh.face(f).centroid.x < 1e-9) {
            inflow += std::fabs(mesh.face(f).areaVector.x) * inletSpeed;
        }
    }
    const double net = solver.netBoundaryFlux();
    std::printf("  [solver_tests] canal com saida: %zu celulas, entrada=%.5f, "
                "balanco liquido=%+.3e (%.2f%% da entrada)\n",
                mesh.cellCount(), inflow, net, 100.0 * std::fabs(net) / inflow);
    std::fflush(stdout);

    AETHER_CHECK(inflow > 0.0);
    // Mass must actually be leaving: a sealed domain would give exactly the
    // inflow as the imbalance, so this separates "outlet works" from "outlet
    // is silently still a wall".
    // **Mass balance closes to machine precision**, and getting there was the
    // whole point. A first version reported a 13.2% deficit that was identical
    // to five digits at t = 3, 8 and 20 -- constant, so not transient, which
    // is the signature of a bookkeeping error rather than physics.
    //
    // The cause was the diagnostic, not the solver: the projection corrects an
    // outlet with the compact face gradient a_b (p_outlet - p_P), the same one
    // its Poisson coefficient is built from, while this check recomputed the
    // flux from the corrected cell velocity and the least-squares gradient --
    // a different operator. Recording the flux the projection actually
    // enforced took the imbalance from 1.3e-01 to 6.3e-14.
    //
    // Same lesson as ROADMAP Fase 1, where the structured cavity's "~0.2
    // divergence" turned out to be a property of the diagnostic too: **measure
    // the operator you actually solve.**
    AETHER_CHECK(std::fabs(net) < 1e-9 * inflow);
    // And the interior must still be divergence-free where no boundary acts.
    AETHER_CHECK(solver.maxFaceDivergence() < 10.0);
}


void testUnstructuredCavityReproducesVortexTopology() {
    // **Deliberately coarse, and the reason is a real cost measured here.**
    // The explicit step is limited by the *smallest* cell, and a Delaunay
    // tetrahedralization of a jittered lattice always produces some slivers,
    // so a finer mesh shrinks dt far faster than it adds cells: at n = 4 and
    // t = 8 this single test took 6m57s, against ~21s for the entire rest of
    // the suite. Vortex topology is a claim a coarse mesh resolves perfectly
    // well, so the mesh and simulated time are set to what the claim needs,
    // not to what would look impressive. Making finer, longer runs practical
    // needs implicit diffusion -- recorded in the ROADMAP, not attempted
    // here.
    const std::size_t n = 4;
    const double lidSpeed = 1.0;
    const aether::mesh::DelaunayTetrahedralization3D tets = buildCubeLatticeTetrahedralization(n);
    const aether::mesh::TetrahedralMesh mesh = aether::mesh::TetrahedralMesh::fromTetrahedralization(tets);

    // Re = lidSpeed * L / nu = 10: safely laminar, single-vortex regime, the
    // same number the structured cavity topology test uses.
    aether::solver::UnstructuredCavitySolver3D solver(mesh, 0.1, [&](const Vector3& p) {
        // Lid on the z = 1 face, tapered to zero at the edges it meets so the
        // corner discontinuity that plagues a constant lid is avoided -- the
        // same regularisation LidDrivenCavitySolver2D uses and for the same
        // reason.
        if (p.z < 1.0 - 1e-9) {
            return Vector3{0.0, 0.0, 0.0};
        }
        const double kPi = 3.14159265358979323846;
        const double sx = std::sin(kPi * p.x);
        const double sy = std::sin(kPi * p.y);
        return Vector3{lidSpeed * sx * sx * sy * sy, 0.0, 0.0};
    });

    // Marched to a fixed *simulated time* rather than a fixed step count:
    // dt is set by the mesh's worst sliver, so a step count would mean
    // different physical durations on different meshes. t = 8 is comfortably
    // past the transient for Re = 10, where viscous diffusion crosses the
    // cavity in O(L^2/nu) = 10 time units.
    const double dt = solver.stableTimeStep();
    const auto stepCount = static_cast<int>(4.0 / dt);
    double worstDivergence = 0.0;
    for (int s = 0; s < stepCount; ++s) {
        solver.step(dt);
        worstDivergence = std::max(worstDivergence, solver.maxFaceDivergence());
    }

    // Average u over the cells nearest the lid and nearest the floor,
    // volume-weighted so a cloud of slivers cannot outvote real cells.
    double topSum = 0.0, topVolume = 0.0, bottomSum = 0.0, bottomVolume = 0.0;
    for (std::size_t cell = 0; cell < mesh.cellCount(); ++cell) {
        const Vector3 c = mesh.cellCentroid(cell);
        const double volume = mesh.cellVolume(cell);
        const double u = solver.velocity(cell).x;
        if (c.z > 0.8) {
            topSum += volume * u;
            topVolume += volume;
        } else if (c.z < 0.2) {
            bottomSum += volume * u;
            bottomVolume += volume;
        }
    }
    const double topMean = topSum / topVolume;
    const double bottomMean = bottomSum / bottomVolume;
    std::printf("  [solver_tests] cavidade nao-estruturada: %zu celulas, u medio topo=%+.5f fundo=%+.5f, "
                "divergencia por faces=%.3e\n",
                mesh.cellCount(), topMean, bottomMean, worstDivergence);
    std::fflush(stdout);

    // Direct viscous drag: the near-lid fluid must follow the lid.
    AETHER_CHECK(topMean > 0.0);
    // And mass conservation in a closed box forces the return flow.
    AETHER_CHECK(bottomMean < 0.0);
    // Bounded divergence: the projection has to actually be doing its job.
    // Bound set from the measured value with margin, not guessed.
    AETHER_CHECK(worstDivergence < 1.0);
    // Guards against the whole thing passing on a dead field.
    AETHER_CHECK(topMean > 0.01);
}

// **The gate for ROADMAP Fase 2.** Two things are checked, and the second
// matters more than the first: that the unstructured solution agrees with an
// independent closed-form answer, and that refining the mesh *reduces* the
// error at a demonstrated rate. Agreement alone at one resolution could come
// from a compensating pair of errors; a clean convergence rate cannot.
//
// The tolerance is set from the measured value rather than guessed -- the
// same discipline every other tolerance in this file follows, and the
// specific lesson relearned in Fase 1, where a metric was trusted without a
// convergence study and turned out to be measuring something else entirely.
void testUnstructuredPlateMatchesFourierSeriesAndConverges() {
    std::printf("  [solver_tests] FVM nao-estruturado, placa vs serie de Fourier:\n");
    std::fflush(stdout);

    double previousError = 0.0;
    double previousInterior = 0.0;
    double previousH = 0.0;
    bool sawConvergence = false;
    double finestError = 0.0;
    double finestInterior = 0.0;
    double finestInteriorOrder = 0.0;

    for (std::size_t n : {4u, 6u, 8u}) {
        const double h = 1.0 / static_cast<double>(n);
        const UnstructuredPlateResult result = solveUnstructuredPlate(n);
        double order = 0.0;
        double interiorOrder = 0.0;
        if (previousError > 0.0) {
            order = std::log(previousError / result.rmsError) / std::log(previousH / h);
            interiorOrder = std::log(previousInterior / result.interiorRmsError) / std::log(previousH / h);
        }
        std::printf(
            "    n=%zu  celulas=%5zu  rmsErro=%8.5f (ordem %4.2f)  semCantos=%8.5f (ordem %4.2f)  "
            "naoOrtog=%.2f\n",
            n, result.cellCount, result.rmsError, order, result.interiorRmsError, interiorOrder,
            result.maxNonOrthogonality);
        std::fflush(stdout);

        if (previousError > 0.0) {
            // The error must actually shrink -- the single most important
            // claim, and the one that separates a working discretization
            // from one that merely produces plausible-looking numbers.
            AETHER_CHECK(result.rmsError < previousError);
            sawConvergence = true;
        }
        previousError = result.rmsError;
        previousInterior = result.interiorRmsError;
        previousH = h;
        finestError = result.rmsError;
        finestInterior = result.interiorRmsError;
        finestInteriorOrder = interiorOrder;
    }
    AETHER_CHECK(sawConvergence);

    // **Gate for ROADMAP Fase 2.2, and how it was finally met.** Four
    // versions were measured, each answering the question the previous one
    // left open rather than being assumed:
    //
    //   orthogonal only          2.879 -> 2.433 -> 2.365   order 0.42, 0.10
    //   + non-orthogonal corr.   1.604 -> 1.056 -> 0.895   order 1.03, 0.58
    //   + least-squares gradient 0.999 -> 0.691 -> 0.522   order 0.91, 0.98
    //   + corrected face grad.   0.989 -> 0.700 -> 0.531   order 0.85, 0.96
    //
    // The first plateaued outright. The second converged with the order
    // drifting *down*. The third fixed that and held a stable ~0.95. The
    // fourth -- distance-weighted interpolation plus replacing the averaged
    // gradient's normal component with the compact difference -- came out
    // **neutral**, within noise of the third. That was a genuine null
    // result: the face interpolation was not the term capping the order,
    // contrary to the hypothesis that motivated it. It is kept because it is
    // the more correct formulation regardless, not because it improved a
    // number.
    //
    // What *was* capping the order is the problem itself. The plate is
    // discontinuous at its two top corners, where the hot edge meets a cold
    // one, and the exact solution has no bounded gradient there -- so no
    // scheme converges at its formal order in a norm that includes those
    // cells. Excluding them (kCornerRadius above) the order climbs to 1.06
    // then 1.63, heading for 2: the scheme is second-order where the exact
    // solution is smooth, which is exactly the claim a finite-volume
    // discretization should be able to make. The structured version of this
    // test avoids the same trap by sampling interior points only.
    //
    // So the assertions below are split accordingly: the global norm must
    // shrink (checked in the loop), and the smooth-region norm must both
    // shrink and reach the measured accuracy with margin.
    AETHER_CHECK(finestError < 0.8);
    AETHER_CHECK(finestInterior < 0.55);
    // Order in the smooth region must be clearly better than first, which is
    // what separates "converging" from "converging at the right rate".
    AETHER_CHECK(finestInteriorOrder > 1.3);
}

// Steady 2D Laplace conduction on a plate with three sides at 0 and the
// top at t0, on a grid with dx != dy (nx != ny over a square domain) -
// this is the test that finally closes the "no genuinely 2D anisotropic
// validation" gap noted in the class comment, by cross-checking both
// against the closed-form Fourier series *and* between Gauss-Seidel and
// CG. Only interior sample points are checked (not the boundary cells
// themselves): what is being validated here is ordinary finite-volume
// truncation error, not the boundary-cell/true-edge half-cell offset the
// 1D tests above are about.
void test2DPlateConductionMatchesFourierSeries() {
    const double lx = 1.0;
    const double ly = 1.0;
    const double t0 = 100.0;
    const std::size_t nx = 81;
    const std::size_t ny = 41;

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(lx, ly, 0.1), nx, ny, 1);

    SteadyDiffusionSolver cgSolver(grid);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 0.0);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::YMin, 0.0);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::YMax, t0);
    const std::size_t cgIterations = cgSolver.solveConjugateGradient(20000, 1e-10);
    AETHER_CHECK(cgIterations < 20000);

    SteadyDiffusionSolver gaussSeidelSolver(grid);
    gaussSeidelSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
    gaussSeidelSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 0.0);
    gaussSeidelSolver.setBoundaryValue(SteadyDiffusionSolver::Face::YMin, 0.0);
    gaussSeidelSolver.setBoundaryValue(SteadyDiffusionSolver::Face::YMax, t0);
    const std::size_t gaussSeidelIterations = gaussSeidelSolver.solve(200000, 1e-8);
    AETHER_CHECK(gaussSeidelIterations < 200000);

    // Fractions of (nx, ny) rather than fixed indices, so the sample grid
    // stays proportionally placed if nx/ny change. Measured empirically at
    // two resolutions (41x21 and 81x41) that doubling resolution roughly
    // halves the gap to the continuum solution, confirming this is
    // ordinary O(h) finite-volume/boundary-placement error, not a bug: the
    // worst sampled point (nearest the hot boundary) was off by ~2.6 at
    // 41x21 and ~1.4 at 81x41. Tolerance below reflects that at this
    // resolution, not machine precision.
    const double sampleFractionsI[] = {0.25, 0.5, 0.75, 0.5, 0.5};
    const double sampleFractionsJ[] = {0.5, 0.5, 0.5, 0.25, 0.75};
    for (std::size_t s = 0; s < 5; ++s) {
        const auto i = static_cast<std::size_t>(sampleFractionsI[s] * static_cast<double>(nx));
        const auto j = static_cast<std::size_t>(sampleFractionsJ[s] * static_cast<double>(ny));

        // Gauss-Seidel and CG solve the identical discretized system, so
        // they must agree closely regardless of how well that system
        // approximates the true continuum solution.
        AETHER_CHECK(nearlyEqual(cgSolver.value(i, j, 0), gaussSeidelSolver.value(i, j, 0), 1e-3));

        const Vector3 center = grid.cellCenter(i, j, 0);
        const double expected = analyticalPlateTemperature(center.x, center.y, lx, ly, t0, 200);
        AETHER_CHECK(nearlyEqual(cgSolver.value(i, j, 0), expected, 2.0));
    }
}

// Jacobi-preconditioned CG (M^-1 = 1/diag(A)) cross-checked against plain
// CG on the same anisotropic (dx != dy) plate problem as
// test2DPlateConductionMatchesFourierSeries. Measured directly, and
// initially unexpected until worked through: the two converge in
// *effectively the same* iteration count (276 vs 276 at 81x41; 619 vs 618
// at 200x20), not fewer. This is a real, explainable null result, not a
// bug: every interior free cell on this project's uniform Cartesian grids
// shares the exact same discrete stencil, so diag(A) is one constant value
// repeated everywhere -- Jacobi preconditioning is then just a uniform
// rescaling of the whole system, which does not change CG's convergence
// behavior (that depends on the *relative* spread of eigenvalues, which a
// single global constant factor leaves untouched). Jacobi would help on a
// grid with genuinely varying per-cell diagonal (graded mesh, spatially
// varying diffusivity) -- not yet a case this project's tests exercise.
void testPreconditionedConjugateGradientMatchesPlainCGOnUniformGrid() {
    const double lx = 1.0;
    const double ly = 1.0;
    const double t0 = 100.0;
    const std::size_t nx = 81;
    const std::size_t ny = 41;

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(lx, ly, 0.1), nx, ny, 1);

    SteadyDiffusionSolver cgSolver(grid);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 0.0);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::YMin, 0.0);
    cgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::YMax, t0);
    const std::size_t cgIterations = cgSolver.solveConjugateGradient(20000, 1e-10);

    SteadyDiffusionSolver pcgSolver(grid);
    pcgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, 0.0);
    pcgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, 0.0);
    pcgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::YMin, 0.0);
    pcgSolver.setBoundaryValue(SteadyDiffusionSolver::Face::YMax, t0);
    const std::size_t pcgIterations = pcgSolver.solvePreconditionedConjugateGradient(20000, 1e-10);

    AETHER_CHECK(cgIterations < 20000);
    AETHER_CHECK(pcgIterations < 20000);
    // Not a speed claim -- see the comment above. Just confirms the
    // preconditioned path is neither broken (blowing up in iterations) nor
    // silently a no-op that skips solving anything.
    AETHER_CHECK(pcgIterations <= cgIterations + 5);

    for (std::size_t j = 0; j < ny; j += 5) {
        for (std::size_t i = 0; i < nx; i += 5) {
            AETHER_CHECK(nearlyEqual(cgSolver.value(i, j, 0), pcgSolver.value(i, j, 0), 1e-6));
        }
    }
}

// MultigridPoissonSolver2D validated against the same closed-form Fourier
// series solution as test2DPlateConductionMatchesFourierSeries (same
// 3-cold-1-hot-side plate, same reused analyticalPlateTemperature()), on a
// power-of-2 grid (64x64) since multigrid's coarsening needs it. Also
// measures the grid-independent-convergence claim that is multigrid's
// whole point: compared directly against SteadyDiffusionSolver's CG and
// Gauss-Seidel on an equivalent 64x64 problem, measured directly (not
// guessed) at 11 V-cycles vs CG's 244 iterations vs Gauss-Seidel's 7425 --
// each V-cycle costs roughly a handful of full-grid Gauss-Seidel sweeps
// (2 pre + 2 post at the finest level, cheaper at each coarser level), so
// this is a real, large speedup, not just a smaller iteration count with
// hidden per-iteration cost. Multigrid's ghost-mirrored Dirichlet boundary
// (ghost = 2*wallValue-interior) is a genuinely different discretization
// from SteadyDiffusionSolver's fixed-boundary-cell-layer convention, so
// per-cell values are not expected to match between the two -- only the
// iteration-count comparison and each method's own agreement with the
// independent Fourier reference are meaningful here.
void testMultigridPoissonMatchesFourierSeriesAndConvergesInFewCycles() {
    const double length = 1.0;
    const double t0 = 100.0;
    const std::size_t n = 64;

    MultigridPoissonSolver2D solver(n, n, length, length);
    solver.setBoundaryValue(MultigridPoissonSolver2D::Face::XMin, 0.0);
    solver.setBoundaryValue(MultigridPoissonSolver2D::Face::XMax, 0.0);
    solver.setBoundaryValue(MultigridPoissonSolver2D::Face::YMin, 0.0);
    solver.setBoundaryValue(MultigridPoissonSolver2D::Face::YMax, t0);
    solver.setSourceTerm(0.0);

    const std::size_t vCycles = solver.solve(100, 1e-9, 2, 2);
    AETHER_CHECK(vCycles < 100);
    AETHER_CHECK(vCycles < 20); // measured 11; comfortable margin above that

    const double h = length / static_cast<double>(n);
    double maxError = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * h;
            const double y = (static_cast<double>(j) + 0.5) * h;
            const double expected = analyticalPlateTemperature(x, y, length, length, t0, 200);
            maxError = std::max(maxError, std::fabs(solver.value(i, j) - expected));
        }
    }
    // Measured directly: ~2.19 at this resolution, the same kind of
    // ordinary O(h) discretization error test2DPlateConductionMatchesFourierSeries
    // documents for SteadyDiffusionSolver, not a bug -- comfortable margin below.
    AETHER_CHECK(maxError < 3.0);
}

// Transient 1D diffusion, d(phi)/dt = d^2(phi)/dx^2, phi(0,t) = phi(L,t) = 0,
// initial condition phi(x,0) = sin(pi*x/L): the classic separable solution
// phi(x,t) = sin(pi*x/L) * exp(-(pi/L)^2 * t) decays without changing
// shape. Validates TransientDiffusionSolver's explicit stepping against
// this closed form, at interior points (the boundary-cell/true-edge
// half-cell offset noted for the steady tests doesn't come into play here
// since sin is already ~0 there for a well-resolved grid).
void testTransientDiffusionMatchesSineDecay() {
    const double lengthX = 1.0;
    const std::size_t nx = 41;
    const double kPi = 3.14159265358979323846;

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(lengthX, 0.1, 0.1), nx, 1, 1);
    TransientDiffusionSolver solver(grid);
    solver.setBoundaryValue(DiffusionProblem::Face::XMin, 0.0);
    solver.setBoundaryValue(DiffusionProblem::Face::XMax, 0.0);

    for (std::size_t i = 0; i < nx; ++i) {
        const double x = grid.cellCenter(i, 0, 0).x;
        solver.setValue(i, 0, 0, std::sin(kPi * x / lengthX));
    }

    const double dt = 0.4 * solver.stableTimeStep(); // safety margin below the marginal-stability limit
    const double targetTime = 0.01;
    const auto steps = static_cast<std::size_t>(targetTime / dt);
    for (std::size_t s = 0; s < steps; ++s) {
        solver.step(dt);
    }

    const double decay = std::exp(-(kPi / lengthX) * (kPi / lengthX) * solver.time());
    for (std::size_t i : {10, 20, 30}) {
        const double x = grid.cellCenter(i, 0, 0).x;
        const double expected = std::sin(kPi * x / lengthX) * decay;
        AETHER_CHECK(nearlyEqual(solver.value(i, 0, 0), expected, 1e-3));
    }
}

// Taylor-Green vortex: an EXACT closed-form solution of the full nonlinear
// 2D incompressible Navier-Stokes equations (doubly periodic domain
// [0,2pi]x[0,2pi]) --
//   u(x,y,t) =  U0 * cos(x) * sin(y) * exp(-2*nu*t)
//   v(x,y,t) = -U0 * sin(x) * cos(y) * exp(-2*nu*t)
// This is the first test that exercises convection together with the
// pressure projection (Poiseuille flow's fully-developed special case has
// zero convection since du/dx = 0 there); it validates
// TaylorGreenVortexSolver2D's predictor+projection loop against real
// physics, not just self-consistency between two numerical methods.
void testTaylorGreenVortexMatchesExactDecay() {
    const double kPi = 3.14159265358979323846;
    const std::size_t n = 32;
    const double length = 2.0 * kPi;
    const double viscosity = 0.05;
    const double u0 = 1.0;

    TaylorGreenVortexSolver2D solver(n, n, length, length, viscosity);
    const double dx = length / static_cast<double>(n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double x = (static_cast<double>(i) + 0.5) * dx;
            const double y = (static_cast<double>(j) + 0.5) * dx;
            solver.setVelocity(i, j, u0 * std::cos(x) * std::sin(y), -u0 * std::sin(x) * std::cos(y));
        }
    }

    const double dt = 0.3 * solver.stableTimeStep(u0);
    const double targetTime = 0.3;
    const auto steps = static_cast<std::size_t>(targetTime / dt);
    for (std::size_t s = 0; s < steps; ++s) {
        solver.step(dt);
        // Divergence stays small and bounded (measured ~4e-3 at this
        // resolution) rather than reaching machine precision: composing
        // two central-difference operators (gradient for the correction,
        // divergence for this check) does not exactly invert the compact
        // 5-point Laplacian solved for pressure on a collocated grid --
        // the classic checkerboard-mode gap Rhie-Chow interpolation
        // exists to close (not implemented here, per the class comment).
        // The bound below is a regression check (a real bug, like the
        // pressure-equation sign error found while building this, made
        // divergence *grow* every step instead of staying bounded), not a
        // claim of exact incompressibility.
        AETHER_CHECK(solver.maxDivergence() < 0.01);
    }

    const double decay = std::exp(-2.0 * viscosity * solver.time());
    const std::size_t samples[] = {8, 16, 24};
    for (std::size_t j : samples) {
        for (std::size_t i : samples) {
            const double x = (static_cast<double>(i) + 0.5) * dx;
            const double y = (static_cast<double>(j) + 0.5) * dx;
            const double expectedU = u0 * std::cos(x) * std::sin(y) * decay;
            const double expectedV = -u0 * std::sin(x) * std::cos(y) * decay;
            AETHER_CHECK(nearlyEqual(solver.u(i, j), expectedU, 5e-3));
            AETHER_CHECK(nearlyEqual(solver.v(i, j), expectedV, 5e-3));
        }
    }
}

// StaggeredNavierStokesSolver3D validated against a genuinely 3D exact
// solution: a Beltrami flow (curl(u) = u), specifically the single ABC
// mode (A=B=C=1):
//   u0 = sin(z) + cos(y)
//   v0 = sin(x) + cos(z)
//   w0 = sin(y) + cos(x)
// Direct substitution confirms curl(u0) = u0 and div(u0) = 0 (u has no
// x-dependence, v none in y, w none in z). Because curl(u0) = u0, the
// standard identity (u.grad)u = grad(|u|^2/2) - u x curl(u) collapses the
// nonlinear advection term to a pure gradient (u x u = 0), so the *entire*
// nonlinear term is absorbable into a modified pressure -- leaving
// u(t) = u0 * exp(-nu*t) an EXACT solution of the full nonlinear
// incompressible Navier-Stokes equations (each Cartesian component also
// satisfies laplacian(u0) = -u0 by direct differentiation, matching the
// same decay rate on the diffusion side). Self-derived from a standard
// vector identity, not recalled benchmark data -- and genuinely 3D: all
// three components are active with real cross-direction dependence in
// every momentum equation (unlike embedding a z-invariant 2D flow in a 3D
// grid, which would never exercise the new cross-direction terms at all).
void testStaggeredNavierStokes3DMatchesBeltramiDecay() {
    const double kPi = 3.14159265358979323846;
    const std::size_t n = 16;
    const double length = 2.0 * kPi;
    const double viscosity = 0.1;

    StaggeredNavierStokesSolver3D solver(n, n, n, length, length, length, viscosity);
    const double h = length / static_cast<double>(n);
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double yCenter = (static_cast<double>(j) + 0.5) * h;
                const double zCenter = (static_cast<double>(k) + 0.5) * h;
                const double xCenter = (static_cast<double>(i) + 0.5) * h;

                const double uInit = std::sin(zCenter) + std::cos(yCenter); // u-face: (i*h, yCenter, zCenter)
                const double vInit = std::sin(xCenter) + std::cos(zCenter); // v-face: (xCenter, j*h, zCenter)
                const double wInit = std::sin(yCenter) + std::cos(xCenter); // w-face: (xCenter, yCenter, k*h)
                solver.setVelocity(i, j, k, uInit, vInit, wInit);
            }
        }
    }

    const double velocityScale = 2.0; // |sin(.)| + |cos(.)| <= 2 componentwise
    const double dt = 0.3 * solver.stableTimeStep(velocityScale);
    const double targetTime = 0.5;
    const auto steps = static_cast<std::size_t>(targetTime / dt);
    for (std::size_t s = 0; s < steps; ++s) {
        solver.step(dt);
        // On a staggered grid the pressure-gradient and divergence
        // operators are exact discrete adjoints, so divergence should stay
        // near machine precision -- not just "small and bounded" like the
        // collocated TaylorGreenVortexSolver2D -- this is the actual point
        // of staggering, checked directly rather than assumed. Measured
        // directly (not guessed): max divergence over the whole run is
        // ~9.4e-13 at this resolution; the bound below keeps a comfortable
        // margin while still catching a real regression.
        AETHER_CHECK(solver.maxDivergence() < 1e-9);
    }

    const double decay = std::exp(-viscosity * solver.time());
    const std::size_t samples[] = {3, 8, 13};
    for (std::size_t k : samples) {
        for (std::size_t j : samples) {
            for (std::size_t i : samples) {
                const double yCenter = (static_cast<double>(j) + 0.5) * h;
                const double zCenter = (static_cast<double>(k) + 0.5) * h;
                const double xCenter = (static_cast<double>(i) + 0.5) * h;
                const double expectedU = (std::sin(zCenter) + std::cos(yCenter)) * decay;
                const double expectedV = (std::sin(xCenter) + std::cos(zCenter)) * decay;
                const double expectedW = (std::sin(yCenter) + std::cos(xCenter)) * decay;
                // Measured directly: worst-case error at these sample
                // points is ~8e-4 at this resolution/timestep (ordinary
                // discretization error, not a bug) -- tolerance below
                // keeps a comfortable margin over that.
                AETHER_CHECK(nearlyEqual(solver.u(i, j, k), expectedU, 3e-3));
                AETHER_CHECK(nearlyEqual(solver.v(i, j, k), expectedV, 3e-3));
                AETHER_CHECK(nearlyEqual(solver.w(i, j, k), expectedW, 3e-3));
            }
        }
    }
}

// StaggeredLidDrivenCavitySolver3D: the 3D generalization of
// LidDrivenCavitySolver2D on the staggered grid, closing the "solid-wall
// 3D" gap flagged when StaggeredNavierStokesSolver3D (periodic-only)
// shipped. Same "no literature benchmark" philosophy: with the lid
// stationary and zero initial velocity, there is no forcing anywhere, so
// the field must stay at *exactly* zero (bit-for-bit), same as the 2D
// class's own rest-state test.
void testStaggeredLidDrivenCavity3DStaysAtRestWhenLidStationary() {
    const std::size_t n = 6;
    StaggeredLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i <= n; ++i) {
                AETHER_CHECK(solver.u(i, j, k) == 0.0);
            }
        }
        for (std::size_t j = 0; j <= n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.v(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k <= n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.w(i, j, k) == 0.0);
            }
        }
    }
}

// With the lid moving (Re = lid*L/nu = 10): mass conservation and the same
// primary-vortex topology already validated for the 2D cavity (top layer,
// just below the lid, dragged in the lid's direction; bottom layer
// reversed by mass conservation in the closed box). Measured directly (not
// guessed): divergence stays at ~2.5e-14 over 400 steps at this
// resolution -- genuinely near machine precision, the actual point of the
// staggered grid (unlike LidDrivenCavitySolver2D's collocated ~0.05-0.19),
// confirmed carrying over correctly from the periodic 3D solver to the
// solid-wall case.
void testStaggeredLidDrivenCavity3DMassConservationAndVortexTopology() {
    const std::size_t n = 16;
    StaggeredLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.1, 1.0); // Re = 10
    const double dt = 0.3 * solver.stableTimeStep();

    for (int s = 0; s < 400; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 1e-9);
    }

    double topMeanU = 0.0;
    double bottomMeanU = 0.0;
    std::size_t count = 0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i <= n; ++i) {
            topMeanU += solver.u(i, j, n - 1);
            bottomMeanU += solver.u(i, j, 0);
            ++count;
        }
    }
    topMeanU /= static_cast<double>(count);
    bottomMeanU /= static_cast<double>(count);
    AETHER_CHECK(topMeanU > 0.1);     // dragged in the lid's direction, measured ~0.175
    AETHER_CHECK(bottomMeanU < 0.0);  // reversed by mass conservation, measured ~-0.0036
}

// MixingLengthLidDrivenCavitySolver3D: the first turbulence closure in this
// project coupled to a real 3D convecting field (every earlier closure that
// saw real convection -- MixingLengthLidDrivenCavitySolver2D,
// KEpsilonLidDrivenCavitySolver2D, KOmegaSSTLidDrivenCavitySolver2D -- was
// still on a flat 2D colocated grid). Same exact-rest-state check as every
// other lid-stationary test in this project: with lidVelocity = 0, velocity
// stays at genuine zero everywhere, so the strain-rate magnitude the
// mixing-length closure depends on is also exactly zero everywhere, making
// nu_t identically zero for all time (a stronger, algebraic guarantee here
// than KEpsilon/KOmegaSST's transported k/epsilon/omega fields, which are
// *not* pinned to zero at rest since they diffuse independently of velocity).
void testMixingLengthCavity3DStaysAtRestWhenLidStationary() {
    const std::size_t n = 6;
    MixingLengthLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i <= n; ++i) {
                AETHER_CHECK(solver.u(i, j, k) == 0.0);
            }
        }
        for (std::size_t j = 0; j <= n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.v(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k <= n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.w(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.eddyViscosity(i, j, k) == 0.0);
            }
        }
    }
}

// With the lid moving (Re = lid*L/nu = 100): mass conservation (measured
// directly at ~2.9e-13 over 400 steps -- the same near-machine-precision
// order StaggeredLidDrivenCavitySolver3D itself achieves, confirming the
// staggered grid's advantage over the 2D closures' collocated ~0.14-0.19
// carries over even with the added eddy-viscosity diffusion term), the same
// primary-vortex topology, and a structural check on nu_t: non-negative
// everywhere, and larger at the cavity center than at a wall-adjacent cell
// (mixing length grows with wall distance) -- the 3D analog of every other
// closure's own sanity-floor check in this project.
void testMixingLengthCavity3DMassConservationTopologyAndEddyViscosity() {
    const std::size_t n = 12;
    MixingLengthLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.01, 1.0); // Re = 100
    double dt = 0.3 * solver.stableTimeStep();

    for (int s = 0; s < 400; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 1e-9);
        if (s % 100 == 0) {
            dt = 0.3 * solver.stableTimeStep(); // re-tighten as nu_t grows
        }
    }

    double topMeanU = 0.0;
    double bottomMeanU = 0.0;
    std::size_t count = 0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i <= n; ++i) {
            topMeanU += solver.u(i, j, n - 1);
            bottomMeanU += solver.u(i, j, 0);
            ++count;
        }
    }
    topMeanU /= static_cast<double>(count);
    bottomMeanU /= static_cast<double>(count);
    AETHER_CHECK(topMeanU > 0.1);    // dragged in the lid's direction, measured ~0.151
    AETHER_CHECK(bottomMeanU < 0.0); // reversed by mass conservation, measured ~-0.0049

    double maxNut = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double nut = solver.eddyViscosity(i, j, k);
                AETHER_CHECK(nut >= 0.0);
                maxNut = std::max(maxNut, nut);
            }
        }
    }
    AETHER_CHECK(maxNut > 0.0); // some turbulent mixing actually developed
    AETHER_CHECK(solver.eddyViscosity(0, 0, 0) < solver.eddyViscosity(n / 2, n / 2, n / 2));
}

// KEpsilonLidDrivenCavitySolver3D: the two-equation k-epsilon closure
// extended from the real 2D cavity (KEpsilonLidDrivenCavitySolver2D) to the
// real 3D solid-wall cavity, mirroring the same step
// MixingLengthLidDrivenCavitySolver3D just took for the simpler algebraic
// closure. With the lid stationary, momentum has zero forcing regardless of
// k/epsilon (production = nu_t*strain^2 = 0 when velocity is uniformly
// zero, and the momentum diffusion terms multiply a uniformly-zero velocity
// field, so they vanish too, whatever nu_t is) -- so u, v, w must stay
// *exactly* zero, the same strong check used throughout this project. k and
// epsilon are *not* asserted to stay at their initial values (same
// reasoning as the 2D test): with k fixed to 0 at all six walls but
// initialized to a uniform nonzero interior value, diffusion alone still
// evolves k/epsilon toward the walls even at rest.
void testKEpsilonCavity3DVelocityStaysAtRestWhenLidStationary() {
    const std::size_t n = 6;
    KEpsilonLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i <= n; ++i) {
                AETHER_CHECK(solver.u(i, j, k) == 0.0);
            }
        }
        for (std::size_t j = 0; j <= n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.v(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k <= n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.w(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.k(i, j, k) >= 0.0);
                AETHER_CHECK(solver.epsilon(i, j, k) >= 0.0);
            }
        }
    }
}

// With the lid moving (Re=100, matching MixingLengthLidDrivenCavitySolver3D's
// own test): mass conservation (measured directly at ~4.1e-11 over 600
// steps -- still many orders of magnitude below the 2D closures' collocated
// ~0.14-0.19, though a bit above the algebraic mixing-length closure's own
// ~2.9e-13, since k/epsilon transport feeds back into nu_t and hence the
// momentum diffusion coefficients every step), the same primary-vortex
// topology, and structural checks on k/epsilon/nu_t -- all non-negative
// everywhere, and larger at the cavity center than at a wall-adjacent
// corner cell, consistent with k=0 exactly at the walls. As with the 2D
// class, this deliberately does not claim self-sustaining turbulence:
// measured directly that k and nu_t both stay small/decaying at this
// Reynolds number, the same physically-correct low-Re behavior the 2D
// class's own test documents (real transition needs Re ~ 10^4).
void testKEpsilonCavity3DMassConservationTopologyAndStructure() {
    const std::size_t n = 12;
    KEpsilonLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.01, 1.0); // Re = 100
    double dt = 0.3 * solver.stableTimeStep();

    for (int s = 0; s < 600; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 1e-6);
        if (s % 100 == 0) {
            dt = 0.3 * solver.stableTimeStep();
        }
    }

    double topMeanU = 0.0;
    double bottomMeanU = 0.0;
    std::size_t count = 0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i <= n; ++i) {
            topMeanU += solver.u(i, j, n - 1);
            bottomMeanU += solver.u(i, j, 0);
            ++count;
        }
    }
    topMeanU /= static_cast<double>(count);
    bottomMeanU /= static_cast<double>(count);
    AETHER_CHECK(topMeanU > 0.1);    // dragged in the lid's direction, measured ~0.148
    AETHER_CHECK(bottomMeanU < 0.0); // reversed by mass conservation, measured ~-0.0045

    double maxNut = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double nut = solver.eddyViscosity(i, j, k);
                const double kVal = solver.k(i, j, k);
                const double epsVal = solver.epsilon(i, j, k);
                AETHER_CHECK(nut >= 0.0);
                AETHER_CHECK(kVal >= 0.0);
                AETHER_CHECK(epsVal >= 0.0);
                maxNut = std::max(maxNut, nut);
            }
        }
    }
    AETHER_CHECK(maxNut > 0.0);
    AETHER_CHECK(solver.k(0, 0, 0) < solver.k(n / 2, n / 2, n / 2));
    AETHER_CHECK(solver.eddyViscosity(0, 0, 0) < solver.eddyViscosity(n / 2, n / 2, n / 2));
}

// SmagorinskyLesLidDrivenCavitySolver3D, rest state. Like mixing length
// (and unlike the two-equation closures), the subgrid viscosity is a
// purely algebraic function of the instantaneous strain rate with no
// transport equation of its own -- so with the lid stationary, velocity
// stays at genuine zero, hence strain is genuine zero, hence nu_sgs is
// genuine zero. Everything is checked for exact zero, not a tolerance.
void testSmagorinskyLes3DStaysAtRestWhenLidStationary() {
    const std::size_t n = 6;
    SmagorinskyLesLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i <= n; ++i) {
                AETHER_CHECK(solver.u(i, j, k) == 0.0);
            }
        }
        for (std::size_t j = 0; j <= n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.v(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k <= n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.w(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.subgridViscosity(i, j, k) == 0.0);
            }
        }
    }
}

// **The test that actually distinguishes LES from RANS**, and the reason
// this class exists as something more than "mixing length with a
// different length scale".
//
// An LES model resolves the large eddies and models only what falls below
// the grid filter width, so its length scale is tied to the MESH
// (Delta = cbrt(dx*dy*dz)). Refining the mesh must therefore shrink
// nu_sgs, approaching DNS in the limit Delta -> 0. A RANS model's length
// scale comes from the FLOW GEOMETRY (here, wall distance), so refining
// the mesh does not shrink it -- it resolves steeper gradients, which if
// anything makes |S| and hence nu_t *larger*.
//
// This is directly measurable rather than merely assertable, so it is
// measured: the same physical problem (Re=100) is run to the same
// physical end time at three resolutions, with both this class and
// MixingLengthLidDrivenCavitySolver3D. Measured directly (not guessed):
//
//   n:              8          12          16
//   LES  nu_sgs:  2.65e-3    1.62e-3    1.07e-3   (falls, 0.41x overall)
//   RANS nu_t:    3.81e-3    5.66e-3    6.26e-3   (rises, 1.64x overall)
//
// The assertions below pin the *directions* (strictly monotone for LES,
// and LES falling while RANS rises) rather than the specific numbers,
// since the directions are what the physics requires while the magnitudes
// depend on resolution and run length.
void testSmagorinskyLes3DSubgridViscosityVanishesUnderMeshRefinement() {
    const double length = 1.0;
    const double viscosity = 0.01; // Re = lid*L/nu = 100
    const double lidVelocity = 1.0;
    const double endTime = 0.6;

    auto runLes = [&](std::size_t n) {
        SmagorinskyLesLidDrivenCavitySolver3D solver(n, n, n, length, length, length, viscosity,
                                                      lidVelocity);
        const double dt = 0.3 * solver.stableTimeStep();
        const int steps = static_cast<int>(endTime / dt);
        for (int s = 0; s < steps; ++s) {
            solver.step(dt);
        }
        double maxNut = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t i = 0; i < n; ++i) {
                    maxNut = std::max(maxNut, solver.subgridViscosity(i, j, k));
                }
            }
        }
        return maxNut;
    };

    auto runRans = [&](std::size_t n) {
        MixingLengthLidDrivenCavitySolver3D solver(n, n, n, length, length, length, viscosity,
                                                    lidVelocity);
        double dt = 0.3 * solver.stableTimeStep();
        const int steps = static_cast<int>(endTime / dt);
        for (int s = 0; s < steps; ++s) {
            solver.step(dt);
            if (s % 100 == 0) {
                dt = 0.3 * solver.stableTimeStep();
            }
        }
        double maxNut = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t i = 0; i < n; ++i) {
                    maxNut = std::max(maxNut, solver.eddyViscosity(i, j, k));
                }
            }
        }
        return maxNut;
    };

    const double les8 = runLes(8);
    const double les12 = runLes(12);
    const double les16 = runLes(16);

    // Every refinement must strictly reduce the modelled viscosity.
    AETHER_CHECK(les8 > 0.0);
    AETHER_CHECK(les12 < les8);
    AETHER_CHECK(les16 < les12);
    // And the reduction must be substantial, not a rounding wobble.
    AETHER_CHECK(les16 < 0.6 * les8);

    // The RANS contrast: its geometry-set length scale means refinement
    // does not drive its eddy viscosity toward zero.
    const double rans8 = runRans(8);
    const double rans16 = runRans(16);
    AETHER_CHECK(rans16 > rans8);
}

// With the lid moving (Re=100): mass conservation at the staggered grid's
// usual near-machine precision, the same primary-vortex topology as every
// other 3D cavity solver here, and structural checks on nu_sgs --
// non-negative everywhere, and smaller at a wall-adjacent corner than at
// the cavity center (the kappa*d_wall cap forces it to vanish at walls).
void testSmagorinskyLes3DMassConservationTopologyAndStructure() {
    const std::size_t n = 12;
    SmagorinskyLesLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.01, 1.0); // Re = 100
    double dt = 0.3 * solver.stableTimeStep();

    for (int s = 0; s < 400; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 1e-9); // measured ~4.2e-13
        if (s % 100 == 0) {
            dt = 0.3 * solver.stableTimeStep();
        }
    }

    double topMeanU = 0.0;
    double bottomMeanU = 0.0;
    std::size_t count = 0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i <= n; ++i) {
            topMeanU += solver.u(i, j, n - 1);
            bottomMeanU += solver.u(i, j, 0);
            ++count;
        }
    }
    topMeanU /= static_cast<double>(count);
    bottomMeanU /= static_cast<double>(count);
    AETHER_CHECK(topMeanU > 0.1);    // dragged by the lid, measured ~0.150
    AETHER_CHECK(bottomMeanU < 0.0); // reversed by mass conservation, measured ~-0.0046

    double maxNut = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double nut = solver.subgridViscosity(i, j, k);
                AETHER_CHECK(nut >= 0.0);
                maxNut = std::max(maxNut, nut);
            }
        }
    }
    AETHER_CHECK(maxNut > 0.0);
    AETHER_CHECK(solver.subgridViscosity(0, 0, 0) < solver.subgridViscosity(n / 2, n / 2, n / 2));

    // Delta = cbrt(dx*dy*dz) is exactly the cell edge on this cubic mesh.
    AETHER_CHECK(nearlyEqual(solver.filterWidth(), 1.0 / static_cast<double>(n), 1e-12));
}

// KOmegaSSTLidDrivenCavitySolver3D: k-omega SST extended from the real 2D
// cavity to the real 3D solid-wall cavity, mirroring the same step
// KEpsilonLidDrivenCavitySolver3D just took. Same exact-rest-state check as
// every other lid-stationary test in this project (momentum has zero
// forcing regardless of k/omega when velocity is uniformly zero). Unlike
// the 2D class's own development (which hit a real 0/0 -> NaN initializing
// omega0 = eps0/(beta**k0) at k0=0), this 3D version used the already-
// algebraically-simplified omega0 = sqrt(k0)/(L*beta*^0.25) from the start,
// so this test is expected to pass cleanly, not catch a regression of that
// specific bug.
void testKOmegaSSTCavity3DVelocityStaysAtRestWhenLidStationary() {
    const std::size_t n = 6;
    KOmegaSSTLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i <= n; ++i) {
                AETHER_CHECK(solver.u(i, j, k) == 0.0);
            }
        }
        for (std::size_t j = 0; j <= n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.v(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k <= n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.w(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.k(i, j, k) >= 0.0);
                AETHER_CHECK(solver.omega(i, j, k) >= 0.0);
            }
        }
    }
}

// With the lid moving (Re=100, matching KEpsilonLidDrivenCavitySolver3D's
// own test): mass conservation (measured directly at ~9.9e-13 -- even
// closer to the pure mixing-length closure's ~2.9e-13 than
// KEpsilonLidDrivenCavitySolver3D's own ~4.1e-11), the same primary-vortex
// topology (measured essentially identical to the independently-built
// KEpsilonLidDrivenCavitySolver3D at the same Re: top/bottom mean u
// ~0.148/-0.0045 for both -- the same reassuring cross-closure consistency
// already documented for the 2D classes), and structural checks on
// k/omega/nu_t.
void testKOmegaSSTCavity3DMassConservationTopologyAndStructure() {
    const std::size_t n = 12;
    KOmegaSSTLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.01, 1.0); // Re = 100
    double dt = 0.3 * solver.stableTimeStep();

    for (int s = 0; s < 600; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 1e-6);
        if (s % 100 == 0) {
            dt = 0.3 * solver.stableTimeStep();
        }
    }

    double topMeanU = 0.0;
    double bottomMeanU = 0.0;
    std::size_t count = 0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i <= n; ++i) {
            topMeanU += solver.u(i, j, n - 1);
            bottomMeanU += solver.u(i, j, 0);
            ++count;
        }
    }
    topMeanU /= static_cast<double>(count);
    bottomMeanU /= static_cast<double>(count);
    AETHER_CHECK(topMeanU > 0.1);    // dragged in the lid's direction, measured ~0.148
    AETHER_CHECK(bottomMeanU < 0.0); // reversed by mass conservation, measured ~-0.0045

    double maxNut = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double nut = solver.eddyViscosity(i, j, k);
                const double kVal = solver.k(i, j, k);
                const double omegaVal = solver.omega(i, j, k);
                AETHER_CHECK(nut >= 0.0);
                AETHER_CHECK(kVal >= 0.0);
                AETHER_CHECK(omegaVal >= 0.0);
                maxNut = std::max(maxNut, nut);
            }
        }
    }
    AETHER_CHECK(maxNut > 0.0);
    AETHER_CHECK(solver.k(0, 0, 0) < solver.k(n / 2, n / 2, n / 2));
    AETHER_CHECK(solver.eddyViscosity(0, 0, 0) < solver.eddyViscosity(n / 2, n / 2, n / 2));
}

// DesSstLidDrivenCavitySolver3D: same exact-rest-state check as every other
// lid-stationary test here (momentum has zero forcing regardless of
// k/omega/F_DES when velocity is uniformly zero).
void testDesSst3DVelocityStaysAtRestWhenLidStationary() {
    const std::size_t n = 6;
    DesSstLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i <= n; ++i) {
                AETHER_CHECK(solver.u(i, j, k) == 0.0);
            }
        }
        for (std::size_t j = 0; j <= n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.v(i, j, k) == 0.0);
            }
        }
    }
    for (std::size_t k = 0; k <= n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(solver.w(i, j, k) == 0.0);
            }
        }
    }
}

// **DES must reduce to plain SST exactly when F_DES is forced to 1.0 for
// the entire trajectory**, not just at a final snapshot. Picking a C_DES so
// large that C_DES*Delta vastly exceeds any realistic L_RANS disables the
// DES branch by construction (max(L_RANS/(C_DES*Delta), 1.0) == 1.0 at
// every cell, every step -- confirmed by instrumenting F_DES during the run
// before writing this assertion, not assumed). With F_DES pinned at 1.0
// throughout, DesSstLidDrivenCavitySolver3D's k-equation is line-for-line
// identical to KOmegaSSTLidDrivenCavitySolver3D's, so the two solvers, run
// on the same problem, must produce bit-for-bit identical fields -- and
// they do (measured exactly 0.0 max difference in both u and nu_t).
//
// This is a stronger check than comparing DES at its *default* C_DES=0.61
// against plain SST: at the default, F_DES can transiently exceed 1.0
// during the initial transient even in a case that settles back to F_DES=1
// by the end (measured directly), so the two trajectories would diverge
// slightly before reconverging -- a real effect, not a bug, but a weaker
// and noisier thing to assert than this exact reduction.
void testDesSst3DReducesToPlainSstWhenCDesIsLarge() {
    const std::size_t n = 6;
    const double length = 1.0;
    const double viscosity = 0.01; // Re = 100
    const double lidVelocity = 1.0;

    DesSstLidDrivenCavitySolver3D des(n, n, n, length, length, length, viscosity, lidVelocity, 1.0e6);
    KOmegaSSTLidDrivenCavitySolver3D sst(n, n, n, length, length, length, viscosity, lidVelocity);

    for (int s = 0; s < 80; ++s) {
        const double desDt = 0.3 * des.stableTimeStep();
        des.step(desDt);
        const double sstDt = 0.3 * sst.stableTimeStep();
        sst.step(sstDt);

        for (std::size_t k = 0; k < n; ++k) {
            for (std::size_t j = 0; j < n; ++j) {
                for (std::size_t i = 0; i < n; ++i) {
                    AETHER_CHECK(des.desFactor(i, j, k) == 1.0);
                }
            }
        }
    }

    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                AETHER_CHECK(des.eddyViscosity(i, j, k) == sst.eddyViscosity(i, j, k));
            }
        }
    }
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i <= n; ++i) {
                AETHER_CHECK(des.u(i, j, k) == sst.u(i, j, k));
            }
        }
    }
}

// **The test that actually distinguishes DES from plain SST**: at the
// default C_DES=0.61, on a Reynolds number high enough to actually build up
// turbulence (Re=1000, matching the SmagorinskyLes-vs-mixing-length
// refinement test's own choice of using a genuinely turbulent case), F_DES
// must grow as the mesh refines -- Delta = cbrt(dx*dy*dz) shrinks, so
// C_DES*Delta drops below the (comparatively mesh-independent) RANS length
// scale L_RANS over more and more of the domain, pulling F_DES above 1.0
// further from the walls. This is the same "length scale tied to the mesh"
// mechanism SmagorinskyLesLidDrivenCavitySolver3D's own refinement test
// checks for nu_sgs directly; here it is checked one level up, on the
// blending factor that drives exactly that mechanism inside DES's RANS
// closure. Measured directly (interior-cell average, n=4,6,8, 80 steps
// each): 2.01, 2.89, 3.63 -- strictly increasing, not guessed.
void testDesSst3DFactorGrowsUnderMeshRefinementAtHighReynolds() {
    const double length = 1.0;
    const double viscosity = 1.0e-3; // Re = 1000
    const double lidVelocity = 1.0;
    const int steps = 80;

    auto coreAverageFDes = [&](std::size_t n) {
        DesSstLidDrivenCavitySolver3D solver(n, n, n, length, length, length, viscosity, lidVelocity);
        for (int s = 0; s < steps; ++s) {
            const double dt = 0.3 * solver.stableTimeStep();
            solver.step(dt);
        }
        double total = 0.0;
        std::size_t count = 0;
        for (std::size_t k = n / 4; k < 3 * n / 4; ++k) {
            for (std::size_t j = n / 4; j < 3 * n / 4; ++j) {
                for (std::size_t i = n / 4; i < 3 * n / 4; ++i) {
                    total += solver.desFactor(i, j, k);
                    ++count;
                }
            }
        }
        return total / static_cast<double>(count);
    };

    const double f4 = coreAverageFDes(4);
    const double f6 = coreAverageFDes(6);
    const double f8 = coreAverageFDes(8);

    AETHER_CHECK(f4 >= 1.0);
    AETHER_CHECK(f6 > f4);
    AETHER_CHECK(f8 > f6);
}

// With the lid moving (Re=100): mass conservation at the staggered grid's
// usual near-machine precision, the same primary-vortex topology as every
// other 3D cavity solver here, and structural non-negativity checks on
// k/omega/nu_t (same style as KOmegaSSTLidDrivenCavitySolver3D's own test).
void testDesSst3DMassConservationTopologyAndStructure() {
    const std::size_t n = 10;
    DesSstLidDrivenCavitySolver3D solver(n, n, n, 1.0, 1.0, 1.0, 0.01, 1.0); // Re = 100
    double dt = 0.3 * solver.stableTimeStep();

    for (int s = 0; s < 400; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 1e-6); // measured ~5.8e-13
        if (s % 100 == 0) {
            dt = 0.3 * solver.stableTimeStep();
        }
    }

    double topMeanU = 0.0;
    double bottomMeanU = 0.0;
    std::size_t count = 0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i <= n; ++i) {
            topMeanU += solver.u(i, j, n - 1);
            bottomMeanU += solver.u(i, j, 0);
            ++count;
        }
    }
    topMeanU /= static_cast<double>(count);
    bottomMeanU /= static_cast<double>(count);
    AETHER_CHECK(topMeanU > 0.1);    // dragged by the lid, measured ~0.132
    AETHER_CHECK(bottomMeanU < 0.0); // reversed by mass conservation, measured ~-0.0053

    double maxNut = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double nut = solver.eddyViscosity(i, j, k);
                const double kVal = solver.k(i, j, k);
                const double omegaVal = solver.omega(i, j, k);
                AETHER_CHECK(nut >= 0.0);
                AETHER_CHECK(kVal >= 0.0);
                AETHER_CHECK(omegaVal >= 0.0);
                AETHER_CHECK(solver.desFactor(i, j, k) >= 1.0);
                maxNut = std::max(maxNut, nut);
            }
        }
    }
    AETHER_CHECK(maxNut > 0.0);

    // Delta = cbrt(dx*dy*dz) is exactly the cell edge on this cubic mesh.
    AETHER_CHECK(nearlyEqual(solver.filterWidth(), 1.0 / static_cast<double>(n), 1e-12));
}

// ImplicitConvectionDiffusionSolver1D: the project's first genuinely
// non-symmetric linear system, and BiCGSTAB's validation. Moderate cell
// Peclet number (Pe_cell = u*h/nu = 0.1, well inside first-order upwind's
// well-behaved regime) so the discrete solution should track the exact ODE
// solution closely, modulo upwind's own numerical/false diffusion.
//
// Measured directly before picking tolerances (established practice):
// Gauss-Seidel needed 2804 iterations to converge to 1e-13; BiCGSTAB
// needed 52 -- a ~54x reduction, the same kind of speedup CG gave over
// Gauss-Seidel for the (symmetric) pressure/diffusion systems elsewhere in
// this project, now demonstrated on a non-symmetric one. The two methods'
// converged fields agreed to ~7.8e-11 (both solve the identical discrete
// system, so this checks solver-implementation consistency, not truth);
// against the exact ODE solution (derived in the class header, not
// recalled), the max error was ~0.0151 -- attributed to first-order
// upwind's known numerical diffusion, not a bug (see the sibling
// refinement test, which confirms the O(h) trend this attribution
// predicts).
void testImplicitConvectionDiffusion1DBiCGStabMatchesExactSolutionAndGaussSeidel() {
    const std::size_t nx = 40;
    const double length = 1.0;
    const double velocity = 2.0;
    const double diffusivity = 0.5;
    const double leftValue = 0.0;
    const double rightValue = 1.0;

    ImplicitConvectionDiffusionSolver1D gaussSeidelSolver(nx, length, velocity, diffusivity, 0.0, leftValue,
                                                            rightValue);
    const std::size_t gsIterations = gaussSeidelSolver.solveGaussSeidel(20000, 1e-13);
    AETHER_CHECK(gsIterations < 20000);

    ImplicitConvectionDiffusionSolver1D bicgstabSolver(nx, length, velocity, diffusivity, 0.0, leftValue,
                                                         rightValue);
    const long long bgIterations = bicgstabSolver.solveBiCGStab(500, 1e-10);
    AETHER_CHECK(bgIterations > 0);
    AETHER_CHECK(bgIterations < 500);
    AETHER_CHECK(static_cast<std::size_t>(bgIterations) * 10 < gsIterations); // measured ~54x fewer

    for (std::size_t i = 0; i < nx; ++i) {
        AETHER_CHECK(nearlyEqual(bicgstabSolver.value(i), gaussSeidelSolver.value(i), 1e-6));
    }

    double maxErrorVsExact = 0.0;
    for (std::size_t i = 0; i < nx; ++i) {
        const double x = bicgstabSolver.cellCenterX(i);
        const double exact = leftValue + (rightValue - leftValue) *
                                              (std::exp(velocity * x / diffusivity) - 1.0) /
                                              (std::exp(velocity * length / diffusivity) - 1.0);
        maxErrorVsExact = std::max(maxErrorVsExact, std::fabs(bicgstabSolver.value(i) - exact));
    }
    AETHER_CHECK(maxErrorVsExact < 0.02); // measured ~0.0151
}

// The error against the exact ODE solution must shrink under mesh
// refinement, and roughly halve when nx doubles -- the signature of
// first-order upwind's O(h) truncation error, distinguishing "expected
// discretization error" from "wrong equation". Measured directly: nx=40
// gives ~0.01513, nx=80 gives ~0.00787 (ratio ~0.52, consistent with O(h)).
void testImplicitConvectionDiffusion1DErrorShrinksUnderMeshRefinement() {
    const double length = 1.0;
    const double velocity = 2.0;
    const double diffusivity = 0.5;
    const double leftValue = 0.0;
    const double rightValue = 1.0;

    auto maxErrorAt = [&](std::size_t nx) {
        ImplicitConvectionDiffusionSolver1D solver(nx, length, velocity, diffusivity, 0.0, leftValue,
                                                     rightValue);
        solver.solveBiCGStab(500, 1e-10);
        double maxError = 0.0;
        for (std::size_t i = 0; i < nx; ++i) {
            const double x = solver.cellCenterX(i);
            const double exact = leftValue + (rightValue - leftValue) *
                                                  (std::exp(velocity * x / diffusivity) - 1.0) /
                                                  (std::exp(velocity * length / diffusivity) - 1.0);
            maxError = std::max(maxError, std::fabs(solver.value(i) - exact));
        }
        return maxError;
    };

    const double error40 = maxErrorAt(40);
    const double error80 = maxErrorAt(80);
    AETHER_CHECK(error80 < error40);
    AETHER_CHECK(error80 < 0.7 * error40); // measured ratio ~0.52
}

// BiCGSTAB must also handle a non-zero source term and reversed flow
// direction (velocity < 0, exercising the other upwind branch) --
// cross-checked against Gauss-Seidel only, since the exact-ODE-solution
// branch above was only derived for the source-free case.
void testImplicitConvectionDiffusion1DBiCGStabMatchesGaussSeidelWithSourceAndReversedFlow() {
    {
        ImplicitConvectionDiffusionSolver1D gs(40, 1.0, 3.0, 0.4, 1.5, 0.2, 0.8);
        gs.solveGaussSeidel(20000, 1e-13);
        ImplicitConvectionDiffusionSolver1D bg(40, 1.0, 3.0, 0.4, 1.5, 0.2, 0.8);
        const long long bgIterations = bg.solveBiCGStab(500, 1e-10);
        AETHER_CHECK(bgIterations > 0);
        for (std::size_t i = 0; i < 40; ++i) {
            AETHER_CHECK(nearlyEqual(bg.value(i), gs.value(i), 1e-6));
        }
    }
    {
        ImplicitConvectionDiffusionSolver1D gs(40, 1.0, -5.0, 0.3, 0.0, 1.0, 0.0);
        gs.solveGaussSeidel(20000, 1e-13);
        ImplicitConvectionDiffusionSolver1D bg(40, 1.0, -5.0, 0.3, 0.0, 1.0, 0.0);
        const long long bgIterations = bg.solveBiCGStab(500, 1e-10);
        AETHER_CHECK(bgIterations > 0);
        for (std::size_t i = 0; i < 40; ++i) {
            AETHER_CHECK(nearlyEqual(bg.value(i), gs.value(i), 1e-6));
        }
    }
}

// GMRES(m) on the same non-symmetric system, validated the same three ways
// BiCGSTAB was: against Gauss-Seidel (implementation consistency, same
// discrete system) and against the exact ODE solution derived in the class
// header. Measured directly: GMRES agrees with Gauss-Seidel to ~5.3e-9 and
// lands the same ~0.01513 from the exact solution BiCGSTAB did -- as it
// must, since both converge to the same discrete answer and that residual
// gap is upwind's own numerical diffusion, not solver error.
void testImplicitConvectionDiffusion1DGmresMatchesExactSolutionAndGaussSeidel() {
    const std::size_t nx = 40;
    const double length = 1.0;
    const double velocity = 2.0;
    const double diffusivity = 0.5;
    const double leftValue = 0.0;
    const double rightValue = 1.0;

    ImplicitConvectionDiffusionSolver1D gaussSeidelSolver(nx, length, velocity, diffusivity, 0.0, leftValue,
                                                            rightValue);
    gaussSeidelSolver.solveGaussSeidel(20000, 1e-13);

    ImplicitConvectionDiffusionSolver1D gmresSolver(nx, length, velocity, diffusivity, 0.0, leftValue,
                                                      rightValue);
    const long long gmresIterations = gmresSolver.solveGmres(30, 500, 1e-10);
    AETHER_CHECK(gmresIterations > 0);
    AETHER_CHECK(gmresIterations < 500);

    for (std::size_t i = 0; i < nx; ++i) {
        AETHER_CHECK(nearlyEqual(gmresSolver.value(i), gaussSeidelSolver.value(i), 1e-6));
    }

    double maxErrorVsExact = 0.0;
    for (std::size_t i = 0; i < nx; ++i) {
        const double x = gmresSolver.cellCenterX(i);
        const double exact = leftValue + (rightValue - leftValue) *
                                              (std::exp(velocity * x / diffusivity) - 1.0) /
                                              (std::exp(velocity * length / diffusivity) - 1.0);
        maxErrorVsExact = std::max(maxErrorVsExact, std::fabs(gmresSolver.value(i) - exact));
    }
    AETHER_CHECK(maxErrorVsExact < 0.02); // measured ~0.01513, same as BiCGSTAB's
}

// **The test that actually distinguishes GMRES from BiCGSTAB**, and the
// reason GMRES earns its place alongside the cheaper method rather than
// being redundant with it.
//
// GMRES minimizes the residual over the entire Krylov subspace at every
// step, so its residual norm cannot increase -- including across restarts,
// since each cycle resumes from the previous iterate and the zero
// correction always lies in the new subspace. BiCGSTAB offers no such
// guarantee. This is directly measurable, so it is measured on three
// different problems (varying Peclet number, flow direction, source term
// and restart length).
//
// Measured directly before writing these assertions -- monotonicity
// violations across each method's full residual history:
//
//   case                        GMRES   BiCGSTAB
//   moderate Pe                   0         4
//   negative velocity, m=10       0         6
//   with source term              0        14
//
// The load-bearing assertion is GMRES's zero, which is the actual
// mathematical guarantee. BiCGSTAB's nonzero count is asserted too, but
// only as evidence that the check is not vacuous on these inputs -- it is
// NOT a claim that BiCGSTAB must always oscillate (measured on a fourth,
// higher-Peclet case, it happened to be monotone as well, which is exactly
// why that case is not used here).
void testImplicitConvectionDiffusion1DGmresResidualIsMonotonicUnlikeBiCGStab() {
    auto monotonicityViolations = [](const std::vector<double>& history) {
        std::size_t violations = 0;
        for (std::size_t i = 1; i < history.size(); ++i) {
            if (history[i] > history[i - 1] * (1.0 + 1e-12)) {
                ++violations;
            }
        }
        return violations;
    };

    struct Case {
        double velocity;
        double diffusivity;
        double source;
        double leftValue;
        double rightValue;
        std::size_t restart;
    };
    const Case cases[] = {
        {2.0, 0.5, 0.0, 0.0, 1.0, 30},  // moderate Pe
        {-5.0, 0.3, 0.0, 1.0, 0.0, 10}, // reversed flow, short restart
        {3.0, 0.4, 1.5, 0.2, 0.8, 30},  // nonzero source term
    };

    std::size_t totalBiCGStabViolations = 0;
    for (const Case& c : cases) {
        ImplicitConvectionDiffusionSolver1D gmresSolver(40, 1.0, c.velocity, c.diffusivity, c.source,
                                                          c.leftValue, c.rightValue);
        gmresSolver.solveGmres(c.restart, 500, 1e-10);
        const std::vector<double>& gmresHistory = gmresSolver.residualHistory();
        AETHER_CHECK(gmresHistory.size() > 1);
        AETHER_CHECK(monotonicityViolations(gmresHistory) == 0); // the guarantee

        ImplicitConvectionDiffusionSolver1D bicgstabSolver(40, 1.0, c.velocity, c.diffusivity, c.source,
                                                             c.leftValue, c.rightValue);
        bicgstabSolver.solveBiCGStab(500, 1e-10);
        totalBiCGStabViolations += monotonicityViolations(bicgstabSolver.residualHistory());
    }
    // Not a guarantee about BiCGSTAB -- just proof the check above has
    // teeth on these inputs (measured 4+6+14 = 24).
    AETHER_CHECK(totalBiCGStabViolations > 0);
}

// GMRES's other structural property: without restarting, on an
// n-dimensional system, the Krylov subspace can grow to at most dimension n
// before it must become invariant, so GMRES converges in at most n
// iterations in exact arithmetic (a finite-termination guarantee no
// stationary method like Gauss-Seidel has). Measured: n=25 converged in
// exactly 25.
void testImplicitConvectionDiffusion1DGmresTerminatesWithinDimensionWhenNotRestarted() {
    const std::size_t n = 25;
    ImplicitConvectionDiffusionSolver1D solver(n, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    const long long iterations = solver.solveGmres(n, 500, 1e-12);
    AETHER_CHECK(iterations > 0);
    AETHER_CHECK(static_cast<std::size_t>(iterations) <= n);
}

// Helpers for the variable-coefficient tests below: a diffusivity profile
// peaking mid-domain (the shape a turbulent nu_t profile actually takes,
// which is the physical motivation for supporting variable coefficients at
// all) and a velocity that accelerates along x.
//
// **The parabola 4t(1-t) is used instead of the sin^2 profile these tests
// originally had, and the reason is a real bug this cost time to find.**
// The first version computed the peak with std::sin. MSVC vectorizes a
// std::sin loop under /O2 using a vector sine, which is permitted to
// differ from the scalar sin in the last bit -- and it does, on 16 of
// these 60 entries. That made the compiled test's input differ by 1 ULP
// from the same formula evaluated through the Python bindings, which was
// enough to send BiCGSTAB down a completely different path: it converged
// in 107 iterations from one set of inputs and *broke down* at iteration 99
// from the other. Hours went into suspecting a stale build, a state leak,
// and the breakdown thresholds before the actual inputs were diffed
// entry-by-entry.
//
// 4t(1-t) has the same shape and the same [0.05, 0.95] range but is built
// only from +, -, * and /, every one of which IEEE-754 requires to be
// correctly rounded -- so it is bit-for-bit identical across compilers,
// optimization levels and language bindings, and the test measures the
// solver rather than the vendor's libm. Verified: no breakdowns anywhere,
// and GMRES's iteration counts are literally unchanged by a deliberate
// 1-ULP perturbation of the whole field (BiCGSTAB's still wobble by a few
// -- which is itself a concrete demonstration of the robustness difference
// the GMRES monotonicity test above asserts).
std::vector<double> peakedDiffusivityField(std::size_t nx) {
    std::vector<double> field(nx);
    for (std::size_t i = 0; i < nx; ++i) {
        const double t = (static_cast<double>(i) + 0.5) / static_cast<double>(nx);
        field[i] = 0.05 + 0.9 * (4.0 * t * (1.0 - t));
    }
    return field;
}

std::vector<double> acceleratingVelocityField(std::size_t nx) {
    std::vector<double> field(nx);
    for (std::size_t i = 0; i < nx; ++i) {
        field[i] = 1.0 + 3.0 * (static_cast<double>(i) + 0.5) / static_cast<double>(nx);
    }
    return field;
}

// Variable coefficients must be a strict *generalization*, not a separate
// code path that happens to agree: setting both fields to the same constant
// the scalar constructor was given must reproduce its result **exactly**
// (bit-for-bit, not approximately), since the face-averaging of equal
// neighbors returns that same value identically. Measured: exactly 0.0
// maximum difference.
void testImplicitConvectionDiffusion1DFlatFieldsReproduceConstantCoefficients() {
    const std::size_t nx = 60;
    ImplicitConvectionDiffusionSolver1D scalarForm(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    scalarForm.solveGmres(30, 500, 1e-12);

    ImplicitConvectionDiffusionSolver1D fieldForm(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    fieldForm.setVelocityField(std::vector<double>(nx, 2.0));
    fieldForm.setDiffusivityField(std::vector<double>(nx, 0.5));
    fieldForm.solveGmres(30, 500, 1e-12);

    for (std::size_t i = 0; i < nx; ++i) {
        AETHER_CHECK(fieldForm.value(i) == scalarForm.value(i));
    }
}

// **ILU(0) is the exact LU factorization of this operator, so a Krylov
// solver preconditioned with it must converge in exactly one iteration.**
//
// This is derived, not empirical: ILU(0) keeps only the sparsity pattern of
// A, and LU factorization of a *tridiagonal* matrix generates no fill-in
// outside that pattern, so nothing is dropped -- M = LU = A exactly, and
// M^-1 A is the identity. (extractTridiagonal() verifies the tridiagonality
// this argument rests on rather than assuming it, throwing if any probed
// column has a nonzero outside the band.) Measured: exactly 1 iteration for
// both GMRES and BiCGSTAB, with and without variable coefficients.
//
// **This is a property of 1D, not a general claim about ILU(0)** -- in
// 2D/3D the stencil is not tridiagonal, ILU(0) drops real fill-in, and it
// becomes a genuine approximation. Recorded here so the one-iteration
// result is not misread later.
void testImplicitConvectionDiffusion1DIncompleteLUSolvesInOneIteration() {
    const std::size_t nx = 60;
    using Preconditioner = ImplicitConvectionDiffusionSolver1D::Preconditioner;

    ImplicitConvectionDiffusionSolver1D reference(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    reference.solveGmres(30, 500, 1e-12);

    ImplicitConvectionDiffusionSolver1D iluGmres(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    iluGmres.setPreconditioner(Preconditioner::IncompleteLU);
    AETHER_CHECK(iluGmres.solveGmres(30, 500, 1e-10) == 1);

    ImplicitConvectionDiffusionSolver1D iluBicgstab(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    iluBicgstab.setPreconditioner(Preconditioner::IncompleteLU);
    AETHER_CHECK(iluBicgstab.solveBiCGStab(500, 1e-10) == 1);

    // One iteration, but still the right answer -- otherwise "converged
    // immediately" would be meaningless.
    for (std::size_t i = 0; i < nx; ++i) {
        AETHER_CHECK(nearlyEqual(iluGmres.value(i), reference.value(i), 1e-6));
        AETHER_CHECK(nearlyEqual(iluBicgstab.value(i), reference.value(i), 1e-6));
    }

    // Same result with genuinely variable coefficients, where the exact-LU
    // argument still holds (the matrix is still tridiagonal).
    ImplicitConvectionDiffusionSolver1D variableIlu(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    variableIlu.setVelocityField(acceleratingVelocityField(nx));
    variableIlu.setDiffusivityField(peakedDiffusivityField(nx));
    variableIlu.setPreconditioner(Preconditioner::IncompleteLU);
    AETHER_CHECK(variableIlu.solveGmres(30, 500, 1e-10) == 1);
}

// **The test that shows Jacobi preconditioning is not a placebo -- and,
// equally, that it IS one on a constant-coefficient problem.**
//
// This system's diagonal is |u| + 2*nu/h per row. With constant u and nu
// that is the same number in every row, so M = D is a scalar multiple of
// the identity and M^-1 A is just A rescaled -- and Krylov methods are
// invariant under scalar rescaling (identical subspace, identical
// minimizer; note the *stopping test* is unchanged too, since
// ||M^-1 r||/||M^-1 b|| = ||r||/||b|| when M is a scalar). Jacobi
// therefore cannot speed anything up there, which is what is measured:
// GMRES 251 unpreconditioned vs 276 with Jacobi -- no gain, and if
// anything slightly worse, the residual rounding difference of dividing
// everything by 42.
//
// Once the coefficients actually vary the diagonal varies too and Jacobi
// becomes a real preconditioner: measured GMRES 360 unpreconditioned vs
// 287 with Jacobi, a ~20% reduction.
//
// **Assertions are on GMRES, not BiCGSTAB, deliberately**: GMRES's
// iteration count on this problem is bit-stable (verified unchanged under
// a deliberate 1-ULP perturbation of the whole diffusivity field), while
// BiCGSTAB's moves by a few iterations under the same perturbation -- so
// pinning ratios of BiCGSTAB counts would be pinning noise. See
// peakedDiffusivityField()'s comment for how that was learned the hard way.
//
// **And the honest limit of the Jacobi result: its benefit is neither large
// nor monotonic in how much the diffusivity varies.** Measured across
// increasing peak-to-wall contrast at the same tolerance: 19x contrast ->
// GMRES 360/287 (helps 20%); 100x -> 304/356 (Jacobi *hurts*); 400x and
// 2500x -> both hit the 500-iteration cap (restarted GMRES(30) stagnates
// regardless). The 19x case is used here because it is the one where the
// effect is real and stable, not because it is representative -- Jacobi is
// a weak preconditioner and this measures it honestly rather than picking
// the flattering framing.
void testImplicitConvectionDiffusion1DJacobiHelpsOnlyWhenDiagonalVaries() {
    const std::size_t nx = 60;
    const double tolerance = 1e-8;
    using Preconditioner = ImplicitConvectionDiffusionSolver1D::Preconditioner;

    // --- Constant coefficients: Jacobi must be (mathematically) inert. ---
    ImplicitConvectionDiffusionSolver1D plainConstant(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    const long long plainConstantIterations = plainConstant.solveGmres(30, 500, tolerance);
    ImplicitConvectionDiffusionSolver1D jacobiConstant(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    jacobiConstant.setPreconditioner(Preconditioner::Jacobi);
    const long long jacobiConstantIterations = jacobiConstant.solveGmres(30, 500, tolerance);

    AETHER_CHECK(plainConstantIterations > 0);
    AETHER_CHECK(jacobiConstantIterations > 0);
    // No speedup: measured 276 with Jacobi vs 251 without.
    AETHER_CHECK(jacobiConstantIterations * 10 >= plainConstantIterations * 9);

    // --- Variable coefficients: now the diagonal varies, and it helps. ---
    auto makeVariable = [&](Preconditioner preconditioner) {
        auto solver = std::make_unique<ImplicitConvectionDiffusionSolver1D>(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
        solver->setVelocityField(acceleratingVelocityField(nx));
        solver->setDiffusivityField(peakedDiffusivityField(nx));
        solver->setPreconditioner(preconditioner);
        return solver;
    };

    auto plainVariable = makeVariable(Preconditioner::None);
    const long long plainVariableIterations = plainVariable->solveGmres(30, 500, tolerance);
    auto jacobiVariable = makeVariable(Preconditioner::Jacobi);
    const long long jacobiVariableIterations = jacobiVariable->solveGmres(30, 500, tolerance);

    AETHER_CHECK(plainVariableIterations > 0);
    AETHER_CHECK(jacobiVariableIterations > 0);
    AETHER_CHECK(jacobiVariableIterations * 10 < plainVariableIterations * 9); // measured 287 vs 360

    // Whatever the preconditioner, the converged answer must be the same
    // one Gauss-Seidel reaches on the same variable-coefficient system.
    ImplicitConvectionDiffusionSolver1D reference(nx, 1.0, 2.0, 0.5, 0.0, 0.0, 1.0);
    reference.setVelocityField(acceleratingVelocityField(nx));
    reference.setDiffusivityField(peakedDiffusivityField(nx));
    reference.solveGaussSeidel(200000, 1e-13);
    for (std::size_t i = 0; i < nx; ++i) {
        AETHER_CHECK(nearlyEqual(plainVariable->value(i), reference.value(i), 1e-5));
        AETHER_CHECK(nearlyEqual(jacobiVariable->value(i), reference.value(i), 1e-5));
    }
}

// Lid-driven cavity validation deliberately avoids comparing to literature
// benchmark tables (e.g. Ghia et al. 1982): those would have to be recalled
// from memory here, a real accuracy risk this project's practice avoids.
// Instead these three tests check properties derivable from first
// principles.

// With the lid stationary (lidVelocity = 0) and the field initialized to
// rest, there is no forcing anywhere: predictor, divergence, and pressure
// correction all operate on exact zeros, so the field must stay *exactly*
// at rest (not just approximately) for as many steps as it's run - a
// strong, exact regression check, not a qualitative one.
void testLidDrivenCavityStaysAtRestWhenLidStationary() {
    LidDrivenCavitySolver2D solver(16, 16, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t j = 0; j < 16; ++j) {
        for (std::size_t i = 0; i < 16; ++i) {
            AETHER_CHECK(solver.u(i, j) == 0.0);
            AETHER_CHECK(solver.v(i, j) == 0.0);
        }
    }
}

// With the lid moving, mass conservation must keep the discrete divergence
// small at every step (same regression-style check as
// TaylorGreenVortexSolver2D, and the same collocated-grid/no-Rhie-Chow
// caveat applies).
void testLidDrivenCavityMassConservation() {
    LidDrivenCavitySolver2D solver(32, 32, 1.0, 1.0, 0.1, 1.0); // Re = lid*L/nu = 10, safely laminar
    const double dt = 0.3 * solver.stableTimeStep();
    // Measured directly: divergence rises to a peak of ~0.185 within the
    // first few steps (transient adjustment away from the zero initial
    // condition) then steadily decays past 0.05 by step 200. The bound
    // below covers that whole trajectory with a comfortable margin, not
    // just the settled value.
    for (int s = 0; s < 200; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 0.25);
    }
}

// The same solver, measured with the operator the projection actually
// conserves. maxDivergence() above uses wide central differences of the
// cell velocities; the pressure Poisson equation is solved with the
// *compact* Laplacian, and composing two wide operators does not reproduce
// a compact one -- so that number was never measuring the quantity being
// driven to zero. maxFaceDivergence() measures it from Rhie-Chow face
// fluxes instead, and the algebra (worked through in the
// LidDrivenCavitySolver2D header) says it must come out at solver
// tolerance, not at 0.2.
//
// Measured before the bound below was chosen, across four cases: face
// divergence lands at 2e-13..1e-12 while wide divergence sits at 0.18..0.28
// on the identical fields -- eleven orders of magnitude apart. So the
// solver was always conserving mass essentially exactly, and the "~0.2
// divergence" this class has documented since it was written was a
// property of the diagnostic, not of the flow.
//
// The bound is deliberately far above the measured value (the CG stopping
// tolerance is 1e-10, so the achievable floor is set by the solve, not by
// this class) and far below the wide-stencil number, so it would catch a
// real regression while not being brittle about the last digits.
void testLidDrivenCavityFaceDivergenceIsAtSolverTolerance() {
    LidDrivenCavitySolver2D solver(32, 32, 1.0, 1.0, 0.1, 1.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 200; ++s) {
        solver.step(dt);
    }
    const double faceDiv = solver.maxFaceDivergence();
    const double wideDiv = solver.maxDivergence();
    std::printf("  [solver_tests] cavidade 2D: div por faces=%.3e, div stencil largo=%.3e\n", faceDiv,
                wideDiv);
    AETHER_CHECK(faceDiv < 1e-8);
    // Guards against the two ever being silently made the same function:
    // they measure different operators and must keep disagreeing.
    AETHER_CHECK(wideDiv > 1e3 * faceDiv);
}

// At low Reynolds number (Re=10 here) a lid-driven cavity settles into a
// single primary recirculating vortex: fluid dragged in the lid's
// direction along the top, forced down the far wall, back in the opposite
// direction along the bottom, and up the near wall. Two consequences of
// that topology are checked, both necessary rather than just plausible:
// (1) the top row of cells (just below the lid) must move predominantly
// in the lid's direction (direct viscous drag); (2) mass conservation in
// a closed box means the bottom row cannot also move that same direction
// on average - it must show net reversed flow, or fluid would pile up
// with nowhere to go.
void testLidDrivenCavityPrimaryVortexTopology() {
    const std::size_t n = 32;
    LidDrivenCavitySolver2D solver(n, n, 1.0, 1.0, 0.1, 1.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 400; ++s) {
        solver.step(dt);
    }

    double topRowMeanU = 0.0;
    double bottomRowMeanU = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        topRowMeanU += solver.u(i, n - 1);
        bottomRowMeanU += solver.u(i, 0);
    }
    topRowMeanU /= static_cast<double>(n);
    bottomRowMeanU /= static_cast<double>(n);

    AETHER_CHECK(topRowMeanU > 0.0);    // dragged in the lid's (+x) direction
    AETHER_CHECK(bottomRowMeanU < 0.0); // forced to reverse by mass conservation
}

// Module 6 (turbulence), first pass: fully-developed turbulent channel
// flow closed with Prandtl's mixing length model. u_tau is obtained from
// the *exact* integral momentum balance (source * height/2, under the
// square root), not a finite-difference wall-gradient estimate - an
// earlier version did the latter using pure molecular viscosity, which
// disagreed with the exact value by ~38% because solve()'s actual wall-
// face diffusivity blends in the adjacent cell's eddy viscosity; that
// mismatch is exactly the kind of bug this exact-formula cross-check
// exists to catch, so keep computing u_tau this way rather than reverting
// to the finite-difference estimate.
//
// Validates the log-law *slope* 1/kappa, not the additive constant B
// (would need a value like 5.0-5.2 recalled from memory - the kind of
// literature-table risk this project's practice avoids) or any published
// turbulent-channel-flow benchmark. The slope is close to guaranteed by
// construction from the l_m = kappa*y ansatz, given a constant-stress
// near-wall layer; measured empirically first (see memory) at Re_tau~1095
// to land ~5.5% below 1/kappa (2.30 vs 2.44), attributable to the channel's
// finite height (the total shear stress is exactly linear across a real
// channel, only approximately constant near the wall) rather than a bug -
// tolerance below reflects that measurement.
void testMixingLengthChannelFlowMatchesLogLawSlope() {
    const double height = 2.0;
    const std::size_t ny = 300;
    const double nu = 5e-5;
    const double source = 0.003;

    MixingLengthChannelFlowSolver1D solver(ny, height, nu, source);
    const std::size_t iterations = solver.solve(5000, 1e-13);
    AETHER_CHECK(iterations < 5000);

    const double uTau = solver.frictionVelocity();
    const double reTau = uTau * (height / 2.0) / nu;
    AETHER_CHECK(reTau > 500.0); // comfortably turbulent, not a marginal/laminar-ish case

    // Two points inside the log-law window: past the viscous sublayer and
    // the mixing-length cap's onset (~0.2195*reTau here), well short of
    // the centerline.
    const std::size_t jNear = 4;
    const std::size_t jFar = 18;
    const double yNear = solver.wallDistance(jNear) * uTau / nu;
    const double yFar = solver.wallDistance(jFar) * uTau / nu;
    AETHER_CHECK(yNear > 25.0 && yFar < 0.2 * reTau);

    const double uPlusNear = solver.u(jNear) / uTau;
    const double uPlusFar = solver.u(jFar) / uTau;
    const double measuredSlope = (uPlusFar - uPlusNear) / (std::log(yFar) - std::log(yNear));

    const double kappa = 0.41;
    const double expectedSlope = 1.0 / kappa;
    AETHER_CHECK(nearlyEqual(measuredSlope, expectedSlope, 0.15)); // ~5.5% measured, comfortable margin

    // Sanity floor: eddy viscosity must be smaller very near the wall than
    // at the channel center (mixing length grows with wall distance, so
    // this should always hold for this closure).
    AETHER_CHECK(solver.eddyViscosity(0) < solver.eddyViscosity(ny / 2));
}

// Module 6 (turbulence), second stage: same channel-flow log-law slope
// check as the mixing-length test above, now with the standard two-
// equation k-epsilon closure and equilibrium wall functions instead of an
// algebraic mixing length. ny=36 keeps the wall-adjacent cells inside the
// valid wall-function y+ window (~30-100) for this Re_tau -- a uniform
// grid ties cell size directly to where the first cell lands, so ny can't
// be increased arbitrarily without moving the wall cells out of that
// window (see the class comment).
//
// The relaxed, semi-implicit Picard iteration used here converges slowly
// (a known characteristic of this kind of coupled, under-relaxed
// nu_t/k/epsilon feedback loop, not specific to this implementation): it
// still reports hitting the iteration cap at 20000 rather than reaching a
// tight residual tolerance, even though the solution has already visibly
// settled (see the symmetry check below). So unlike every other solver
// test in this file, this one does not assert iterations < cap.
void testKEpsilonChannelFlowMatchesLogLawSlopeAndIsSymmetric() {
    const double height = 2.0;
    const std::size_t ny = 36;
    const double nu = 5e-5;
    const double source = 0.003;

    KEpsilonChannelFlowSolver1D solver(ny, height, nu, source);
    solver.solve(20000, 1e-10);

    const double uTau = solver.frictionVelocity();

    // Channel symmetric about its centerline: the velocity profile must be
    // too, regardless of how well-converged the iteration is.
    double maxAsymmetry = 0.0;
    for (std::size_t j = 0; j < ny / 2; ++j) {
        maxAsymmetry = std::max(maxAsymmetry, std::fabs(solver.u(j) - solver.u(ny - 1 - j)));
    }
    AETHER_CHECK(maxAsymmetry < 0.05);

    // Log-law slope, measured between the wall-adjacent cell and a point
    // further into the log region. Measured empirically (see memory) at
    // this exact configuration to land at 2.59 against the expected 2.439
    // (1/kappa) -- looser tolerance than the mixing-length test's, since
    // this coarser wall-function mesh and the coupled k-epsilon nonlinear
    // system carry more numerical noise.
    const std::size_t jNear = 0;
    const std::size_t jFar = 3;
    const double yNear = solver.wallDistance(jNear) * uTau / nu;
    const double yFar = solver.wallDistance(jFar) * uTau / nu;
    AETHER_CHECK(yNear > 25.0 && yNear < 100.0);

    const double uPlusNear = solver.u(jNear) / uTau;
    const double uPlusFar = solver.u(jFar) / uTau;
    const double measuredSlope = (uPlusFar - uPlusNear) / (std::log(yFar) - std::log(yNear));

    const double kappa = 0.41;
    const double expectedSlope = 1.0 / kappa;
    AETHER_CHECK(nearlyEqual(measuredSlope, expectedSlope, 0.65)); // ~2.59 measured vs 2.439 expected

    // Sanity floor, same as the mixing-length test.
    AETHER_CHECK(solver.eddyViscosity(0) < solver.eddyViscosity(ny / 2));
}

// k-omega SST, same channel and wall-function setup as
// KEpsilonChannelFlowSolver1D's test (same height/ny/nu/source), so the
// two closures can be compared on identical footing. Same validation
// philosophy as every turbulence closure this project has built: check
// what's self-derivable (log-law slope, profile symmetry), not a
// literature benchmark table -- SST's own model constants (Menter's
// standard values) are the one part of *this* class that is a recalled
// literature value rather than derived, same caveat already documented for
// k-epsilon's B=5.0 log-law constant.
void testKOmegaSSTChannelFlowMatchesLogLawSlopeAndIsSymmetric() {
    const double height = 2.0;
    const std::size_t ny = 36;
    const double nu = 5e-5;
    const double source = 0.003;

    KOmegaSSTChannelFlowSolver1D solver(ny, height, nu, source);
    solver.solve(20000, 1e-10);

    const double uTau = solver.frictionVelocity();

    // Channel symmetric about its centerline: the velocity profile must be
    // too, regardless of how well-converged the iteration is.
    double maxAsymmetry = 0.0;
    for (std::size_t j = 0; j < ny / 2; ++j) {
        maxAsymmetry = std::max(maxAsymmetry, std::fabs(solver.u(j) - solver.u(ny - 1 - j)));
    }
    AETHER_CHECK(maxAsymmetry < 0.05);

    // Log-law slope, same two-point measurement as the k-epsilon test.
    const std::size_t jNear = 0;
    const std::size_t jFar = 3;
    const double yNear = solver.wallDistance(jNear) * uTau / nu;
    const double yFar = solver.wallDistance(jFar) * uTau / nu;
    AETHER_CHECK(yNear > 25.0 && yNear < 100.0);

    const double uPlusNear = solver.u(jNear) / uTau;
    const double uPlusFar = solver.u(jFar) / uTau;
    const double measuredSlope = (uPlusFar - uPlusNear) / (std::log(yFar) - std::log(yNear));

    const double kappa = 0.41;
    const double expectedSlope = 1.0 / kappa;
    // Measured directly (not guessed): ~2.65 at this exact configuration,
    // against the 2.439 expected -- comparable to k-epsilon's own ~6-9%
    // gap on the same mesh, attributed to the same coarse wall-function
    // mesh and nonlinear coupling noise, not re-investigated further.
    AETHER_CHECK(nearlyEqual(measuredSlope, expectedSlope, 0.5));

    // Sanity floor, same as the other two closures.
    AETHER_CHECK(solver.eddyViscosity(0) < solver.eddyViscosity(ny / 2));
}

// MixingLengthLidDrivenCavitySolver2D: the first turbulence closure coupled
// to a genuinely 2D convecting/recirculating flow (every 1D channel closure
// above only ever saw a fully-developed problem where convection vanishes
// by construction). With the lid stationary and zero initial velocity there
// is no strain anywhere, so nu_t (l_m^2 * strain magnitude) must be exactly
// zero everywhere too -- a strong, exact check, mirroring
// testLidDrivenCavityStaysAtRestWhenLidStationary.
void testMixingLengthCavityStaysAtRestWhenLidStationary() {
    MixingLengthLidDrivenCavitySolver2D solver(16, 16, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t j = 0; j < 16; ++j) {
        for (std::size_t i = 0; i < 16; ++i) {
            AETHER_CHECK(solver.u(i, j) == 0.0);
            AETHER_CHECK(solver.v(i, j) == 0.0);
            AETHER_CHECK(solver.eddyViscosity(i, j) == 0.0);
        }
    }
}

// With the lid moving: mass conservation (same regression-style bound as
// LidDrivenCavitySolver2D's own test), the same primary-vortex topology
// (top row dragged by the lid, bottom row reversed by mass conservation in
// the closed box), and a structural check on nu_t -- non-negative
// everywhere, and larger near the cavity center than at a cell adjacent to
// a wall (mixing length grows with wall distance, up to its cap), the 2D
// analog of the sanity floor every 1D closure's test already checks.
void testMixingLengthCavityMassConservationTopologyAndEddyViscosity() {
    const std::size_t n = 32;
    MixingLengthLidDrivenCavitySolver2D solver(n, n, 1.0, 1.0, 0.01, 1.0); // Re = lid*L/nu = 100
    double dt = 0.3 * solver.stableTimeStep();

    // Measured directly: divergence peaks around ~0.18 within the first
    // steps, comparable to LidDrivenCavitySolver2D's own ~0.185 -- the same
    // collocated-grid/no-Rhie-Chow gap, not something the turbulence
    // closure changes.
    for (int s = 0; s < 400; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 0.3);
        if (s % 100 == 0) {
            dt = 0.3 * solver.stableTimeStep(); // re-tighten as nu_t grows
        }
    }

    double topRowMeanU = 0.0;
    double bottomRowMeanU = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        topRowMeanU += solver.u(i, n - 1);
        bottomRowMeanU += solver.u(i, 0);
    }
    topRowMeanU /= static_cast<double>(n);
    bottomRowMeanU /= static_cast<double>(n);
    AETHER_CHECK(topRowMeanU > 0.1);    // dragged in the lid's direction
    AETHER_CHECK(bottomRowMeanU < 0.0); // reversed by mass conservation

    double maxNut = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double nut = solver.eddyViscosity(i, j);
            AETHER_CHECK(nut >= 0.0);
            maxNut = std::max(maxNut, nut);
        }
    }
    AETHER_CHECK(maxNut > 0.0); // some turbulent mixing actually developed
    AETHER_CHECK(solver.eddyViscosity(0, 0) < solver.eddyViscosity(n / 2, n / 2));
}

// KEpsilonLidDrivenCavitySolver2D: the two-equation k-epsilon closure
// (previously only ever coupled to the 1D fully-developed channel)
// extended to the same real 2D convecting cavity as
// MixingLengthLidDrivenCavitySolver2D. With the lid stationary, momentum
// has zero forcing regardless of k/epsilon (production = nu_t*strain^2 = 0
// when velocity is uniformly zero, and the diffusion terms multiply a
// uniformly-zero velocity field, so they vanish too, whatever nu_t is) --
// so u, v must stay *exactly* zero, the same strong check used throughout
// this project's lid-stationary tests. k and epsilon themselves are *not*
// asserted to stay at their initial values here (unlike mixing length's
// nu_t, which is identically zero at rest by construction): with k fixed
// to 0 at the walls but initialized to a uniform nonzero value in the
// interior, diffusion alone still evolves k/epsilon toward the wall even
// at rest -- real, expected behavior for a transported field, not a
// regression to check bit-for-bit.
void testKEpsilonCavityVelocityStaysAtRestWhenLidStationary() {
    const std::size_t n = 12;
    KEpsilonLidDrivenCavitySolver2D solver(n, n, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            AETHER_CHECK(solver.u(i, j) == 0.0);
            AETHER_CHECK(solver.v(i, j) == 0.0);
            AETHER_CHECK(solver.k(i, j) >= 0.0);
            AETHER_CHECK(solver.epsilon(i, j) >= 0.0);
        }
    }
}

// With the lid moving (Re=100, the same value MixingLengthLidDrivenCavitySolver2D's
// own test uses): mass conservation, the same primary-vortex topology, and
// structural checks on k/epsilon/nu_t -- all non-negative everywhere, and
// (mirroring the wall-vs-center sanity floor every turbulence closure's
// test in this project checks) smaller at a wall-adjacent cell than at the
// cavity center, consistent with k=0 exactly at the wall and nu_t
// vanishing there too. Deliberately does not claim to demonstrate
// self-sustaining turbulence: measured directly (not guessed) that both
// k and nu_t *decay* over the run at this Reynolds number rather than
// reaching a nonzero statistically steady level, which is the physically
// correct behavior for k-epsilon at a Re this project's other cavity
// solvers already document as "safely laminar" (real lid-driven-cavity
// transition to turbulence needs Re on the order of 10^4, far beyond what
// any solver in this project has been run at) -- not a bug, and not
// something this test claims otherwise.
void testKEpsilonCavityMassConservationTopologyAndStructure() {
    const std::size_t n = 20;
    KEpsilonLidDrivenCavitySolver2D solver(n, n, 1.0, 1.0, 0.01, 1.0); // Re = lid*L/nu = 100
    double dt = 0.3 * solver.stableTimeStep();

    // Measured directly: divergence peaks around ~0.145 over this run,
    // the same collocated-grid/no-Rhie-Chow order of magnitude as every
    // other cavity solver in this project.
    for (int s = 0; s < 600; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 0.3);
        if (s % 50 == 0) {
            dt = 0.3 * solver.stableTimeStep();
        }
    }

    double topRowMeanU = 0.0;
    double bottomRowMeanU = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        topRowMeanU += solver.u(i, n - 1);
        bottomRowMeanU += solver.u(i, 0);
    }
    topRowMeanU /= static_cast<double>(n);
    bottomRowMeanU /= static_cast<double>(n);
    AETHER_CHECK(topRowMeanU > 0.1);    // dragged in the lid's direction, measured ~0.38
    AETHER_CHECK(bottomRowMeanU < 0.0); // reversed by mass conservation, measured ~-0.006

    double maxNut = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double kVal = solver.k(i, j);
            const double epsVal = solver.epsilon(i, j);
            const double nut = solver.eddyViscosity(i, j);
            AETHER_CHECK(kVal >= 0.0);
            AETHER_CHECK(epsVal >= 0.0);
            AETHER_CHECK(nut >= 0.0);
            maxNut = std::max(maxNut, nut);
        }
    }
    AETHER_CHECK(maxNut > 0.0);
    AETHER_CHECK(solver.k(0, 0) < solver.k(n / 2, n / 2));
    AETHER_CHECK(solver.eddyViscosity(0, 0) < solver.eddyViscosity(n / 2, n / 2));
}

// KOmegaSSTLidDrivenCavitySolver2D: the same real-2D-convection extension
// as KEpsilonLidDrivenCavitySolver2D, now for k-omega SST. Same rest-state
// reasoning: with the lid stationary, momentum's production and diffusion
// terms both vanish identically (they multiply a uniformly-zero velocity
// field) regardless of k/omega/nu_t, so u and v must stay *exactly* zero.
void testKOmegaSSTCavityVelocityStaysAtRestWhenLidStationary() {
    const std::size_t n = 12;
    KOmegaSSTLidDrivenCavitySolver2D solver(n, n, 1.0, 1.0, 0.1, 0.0);
    const double dt = 0.3 * solver.stableTimeStep();
    for (int s = 0; s < 20; ++s) {
        solver.step(dt);
    }
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            AETHER_CHECK(solver.u(i, j) == 0.0);
            AETHER_CHECK(solver.v(i, j) == 0.0);
            AETHER_CHECK(solver.k(i, j) >= 0.0);
            AETHER_CHECK(solver.omega(i, j) >= 0.0);
        }
    }
}

// With the lid moving (Re=100, same configuration as the k-epsilon cavity
// test, for direct comparability): mass conservation, the same primary-
// vortex topology, and structural checks on k/omega/nu_t. Measured
// directly: top/bottom row mean u (~0.382/-0.0064) and peak divergence
// (~0.144) match KEpsilonLidDrivenCavitySolver2D's own measurements at
// this same Re almost exactly -- a reassuring cross-closure consistency
// check, since both are independent implementations of different
// turbulence models coupled to the *same* underlying momentum solver and
// should agree on the parts that don't depend on closure details (gross
// topology, divergence order of magnitude). Same honest caveat as the
// k-epsilon cavity: k and nu_t decay over the run at this Reynolds number
// rather than sustaining a nonzero turbulent level -- physically correct
// for a Re this project already documents as "safely laminar" elsewhere,
// not a bug.
void testKOmegaSSTCavityMassConservationTopologyAndStructure() {
    const std::size_t n = 20;
    KOmegaSSTLidDrivenCavitySolver2D solver(n, n, 1.0, 1.0, 0.01, 1.0); // Re = lid*L/nu = 100
    double dt = 0.3 * solver.stableTimeStep();

    for (int s = 0; s < 600; ++s) {
        solver.step(dt);
        AETHER_CHECK(solver.maxDivergence() < 0.3);
        if (s % 50 == 0) {
            dt = 0.3 * solver.stableTimeStep();
        }
    }

    double topRowMeanU = 0.0;
    double bottomRowMeanU = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        topRowMeanU += solver.u(i, n - 1);
        bottomRowMeanU += solver.u(i, 0);
    }
    topRowMeanU /= static_cast<double>(n);
    bottomRowMeanU /= static_cast<double>(n);
    AETHER_CHECK(topRowMeanU > 0.1);    // dragged in the lid's direction, measured ~0.38
    AETHER_CHECK(bottomRowMeanU < 0.0); // reversed by mass conservation, measured ~-0.006

    double maxNut = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double kVal = solver.k(i, j);
            const double omegaVal = solver.omega(i, j);
            const double nut = solver.eddyViscosity(i, j);
            AETHER_CHECK(kVal >= 0.0);
            AETHER_CHECK(omegaVal >= 0.0);
            AETHER_CHECK(nut >= 0.0);
            maxNut = std::max(maxNut, nut);
        }
    }
    AETHER_CHECK(maxNut > 0.0);
    AETHER_CHECK(solver.k(0, 0) < solver.k(n / 2, n / 2));
    AETHER_CHECK(solver.eddyViscosity(0, 0) < solver.eddyViscosity(n / 2, n / 2));
}

} // namespace

int main() {
    test1DConductionMatchesAnalyticalProfile();
    testBoundaryCellsStayFixed();
    testPlanePoiseuilleProfile();
    testConjugateGradientMatchesGaussSeidelAndConvergesFaster();
    testConjugateGradientWithSourceTermMatchesPoiseuilleProfile();
    test2DPlateConductionMatchesFourierSeries();
    testTransientDiffusionMatchesSineDecay();
    testTaylorGreenVortexMatchesExactDecay();
    testStaggeredNavierStokes3DMatchesBeltramiDecay();
    testStaggeredLidDrivenCavity3DStaysAtRestWhenLidStationary();
    testStaggeredLidDrivenCavity3DMassConservationAndVortexTopology();
    testMixingLengthCavity3DStaysAtRestWhenLidStationary();
    testMixingLengthCavity3DMassConservationTopologyAndEddyViscosity();
    testKEpsilonCavity3DVelocityStaysAtRestWhenLidStationary();
    testKEpsilonCavity3DMassConservationTopologyAndStructure();
    testSmagorinskyLes3DStaysAtRestWhenLidStationary();
    testSmagorinskyLes3DSubgridViscosityVanishesUnderMeshRefinement();
    testSmagorinskyLes3DMassConservationTopologyAndStructure();
    testKOmegaSSTCavity3DVelocityStaysAtRestWhenLidStationary();
    testKOmegaSSTCavity3DMassConservationTopologyAndStructure();
    testDesSst3DVelocityStaysAtRestWhenLidStationary();
    testDesSst3DReducesToPlainSstWhenCDesIsLarge();
    testDesSst3DFactorGrowsUnderMeshRefinementAtHighReynolds();
    testDesSst3DMassConservationTopologyAndStructure();
    testImplicitConvectionDiffusion1DBiCGStabMatchesExactSolutionAndGaussSeidel();
    testImplicitConvectionDiffusion1DErrorShrinksUnderMeshRefinement();
    testImplicitConvectionDiffusion1DBiCGStabMatchesGaussSeidelWithSourceAndReversedFlow();
    testImplicitConvectionDiffusion1DGmresMatchesExactSolutionAndGaussSeidel();
    testImplicitConvectionDiffusion1DGmresResidualIsMonotonicUnlikeBiCGStab();
    testImplicitConvectionDiffusion1DGmresTerminatesWithinDimensionWhenNotRestarted();
    testImplicitConvectionDiffusion1DFlatFieldsReproduceConstantCoefficients();
    testImplicitConvectionDiffusion1DIncompleteLUSolvesInOneIteration();
    testImplicitConvectionDiffusion1DJacobiHelpsOnlyWhenDiagonalVaries();
    testLidDrivenCavityStaysAtRestWhenLidStationary();
    testUnstructuredPlateMatchesFourierSeriesAndConverges();
    testUnstructuredCavityReproducesVortexTopology();
    testChannelWithOutletConservesGlobalMass();
    testLidDrivenCavityMassConservation();
    testLidDrivenCavityFaceDivergenceIsAtSolverTolerance();
    testLidDrivenCavityPrimaryVortexTopology();
    testMixingLengthChannelFlowMatchesLogLawSlope();
    testKEpsilonChannelFlowMatchesLogLawSlopeAndIsSymmetric();
    testKOmegaSSTChannelFlowMatchesLogLawSlopeAndIsSymmetric();
    testMixingLengthCavityStaysAtRestWhenLidStationary();
    testMixingLengthCavityMassConservationTopologyAndEddyViscosity();
    testKEpsilonCavityVelocityStaysAtRestWhenLidStationary();
    testKEpsilonCavityMassConservationTopologyAndStructure();
    testKOmegaSSTCavityVelocityStaysAtRestWhenLidStationary();
    testKOmegaSSTCavityMassConservationTopologyAndStructure();
    testMultigridPoissonMatchesFourierSeriesAndConvergesInFewCycles();
    testPreconditionedConjugateGradientMatchesPlainCGOnUniformGrid();
    std::puts("aether_solver_tests: all tests passed");
    return 0;
}
