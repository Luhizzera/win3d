#include "aether/solver/DesSstLidDrivenCavitySolver3D.hpp"
#include "aether/solver/DiffusionProblem.hpp"
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
#include "aether/solver/MultigridPoissonSolver2D.hpp"
#include "aether/solver/SmagorinskyLesLidDrivenCavitySolver3D.hpp"
#include "aether/solver/StaggeredLidDrivenCavitySolver3D.hpp"
#include "aether/solver/StaggeredNavierStokesSolver3D.hpp"
#include "aether/solver/SteadyDiffusionSolver.hpp"
#include "aether/solver/TaylorGreenVortexSolver2D.hpp"
#include "aether/solver/TransientDiffusionSolver.hpp"
#include "aether/solver/UnstructuredCavitySolver3D.hpp"
#include "aether/solver/UnstructuredDiffusionSolver.hpp"
#include "aether/solver/UnstructuredScalarTransportSolver.hpp"

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>

namespace py = pybind11;
using namespace aether::solver;

PYBIND11_MODULE(aether_solver_py, m) {
    m.doc() = "Python bindings for the Aether CFD engine's first FVM solver (Module 5)";

    py::enum_<DiffusionProblem::Face>(m, "Face")
        .value("X_MIN", DiffusionProblem::Face::XMin)
        .value("X_MAX", DiffusionProblem::Face::XMax)
        .value("Y_MIN", DiffusionProblem::Face::YMin)
        .value("Y_MAX", DiffusionProblem::Face::YMax)
        .value("Z_MIN", DiffusionProblem::Face::ZMin)
        .value("Z_MAX", DiffusionProblem::Face::ZMax);

    py::class_<DiffusionProblem>(m, "DiffusionProblem")
        .def("set_boundary_value", &DiffusionProblem::setBoundaryValue, py::arg("face"), py::arg("value"))
        .def("set_source_term", &DiffusionProblem::setSourceTerm, py::arg("source"))
        .def("value", &DiffusionProblem::value);

    py::class_<SteadyDiffusionSolver, DiffusionProblem>(m, "SteadyDiffusionSolver")
        .def(py::init<const aether::mesh::StructuredGrid3D&>(), py::arg("grid"),
             py::keep_alive<1, 2>()) // SteadyDiffusionSolver only stores a pointer to the grid
        .def("solve", &SteadyDiffusionSolver::solve, py::arg("max_iterations") = 10000,
             py::arg("tolerance") = 1e-9)
        .def("solve_conjugate_gradient", &SteadyDiffusionSolver::solveConjugateGradient,
             py::arg("max_iterations") = 10000, py::arg("tolerance") = 1e-9)
        .def("solve_preconditioned_conjugate_gradient",
             &SteadyDiffusionSolver::solvePreconditionedConjugateGradient, py::arg("max_iterations") = 10000,
             py::arg("tolerance") = 1e-9);

    py::class_<TransientDiffusionSolver, DiffusionProblem>(m, "TransientDiffusionSolver")
        .def(py::init<const aether::mesh::StructuredGrid3D&>(), py::arg("grid"),
             py::keep_alive<1, 2>()) // TransientDiffusionSolver only stores a pointer to the grid
        .def("set_value", &TransientDiffusionSolver::setValue, py::arg("i"), py::arg("j"), py::arg("k"),
             py::arg("value"))
        .def("stable_time_step", &TransientDiffusionSolver::stableTimeStep)
        .def("step", &TransientDiffusionSolver::step, py::arg("dt"))
        .def("time", &TransientDiffusionSolver::time);

    py::class_<TaylorGreenVortexSolver2D>(m, "TaylorGreenVortexSolver2D")
        .def(py::init<std::size_t, std::size_t, double, double, double>(), py::arg("nx"), py::arg("ny"),
             py::arg("length_x"), py::arg("length_y"), py::arg("viscosity"))
        .def("set_velocity", &TaylorGreenVortexSolver2D::setVelocity, py::arg("i"), py::arg("j"),
             py::arg("u"), py::arg("v"))
        .def("stable_time_step", &TaylorGreenVortexSolver2D::stableTimeStep, py::arg("velocity_scale"))
        .def("step", &TaylorGreenVortexSolver2D::step, py::arg("dt"))
        .def("u", &TaylorGreenVortexSolver2D::u)
        .def("v", &TaylorGreenVortexSolver2D::v)
        .def("pressure", &TaylorGreenVortexSolver2D::pressure)
        .def("time", &TaylorGreenVortexSolver2D::time)
        .def("max_divergence", &TaylorGreenVortexSolver2D::maxDivergence);

    py::class_<LidDrivenCavitySolver2D> lidCavity2D(m, "LidDrivenCavitySolver2D");
    py::enum_<LidDrivenCavitySolver2D::ConvectionScheme>(lidCavity2D, "ConvectionScheme")
        .value("CENTRAL", LidDrivenCavitySolver2D::ConvectionScheme::Central)
        .value("FIRST_ORDER_UPWIND", LidDrivenCavitySolver2D::ConvectionScheme::FirstOrderUpwind)
        .value("LIMITED_LINEAR_UPWIND",
               LidDrivenCavitySolver2D::ConvectionScheme::LimitedLinearUpwind);
    lidCavity2D
        .def(py::init<std::size_t, std::size_t, double, double, double, double,
                      LidDrivenCavitySolver2D::ConvectionScheme>(), py::arg("nx"),
             py::arg("ny"), py::arg("length_x"), py::arg("length_y"), py::arg("viscosity"),
             py::arg("lid_velocity"),
             py::arg("convection") =
                 LidDrivenCavitySolver2D::ConvectionScheme::LimitedLinearUpwind)
        .def("stable_time_step", &LidDrivenCavitySolver2D::stableTimeStep)
        .def("step", &LidDrivenCavitySolver2D::step, py::arg("dt"))
        .def("u", &LidDrivenCavitySolver2D::u)
        .def("v", &LidDrivenCavitySolver2D::v)
        .def("pressure", &LidDrivenCavitySolver2D::pressure)
        .def("time", &LidDrivenCavitySolver2D::time)
        .def("max_divergence", &LidDrivenCavitySolver2D::maxDivergence)
        .def("max_face_divergence", &LidDrivenCavitySolver2D::maxFaceDivergence)
        .def("load_state", &LidDrivenCavitySolver2D::loadState, py::arg("u"), py::arg("v"), py::arg("p"),
             py::arg("time"));

    py::class_<StaggeredNavierStokesSolver3D>(m, "StaggeredNavierStokesSolver3D")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, double, double, double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("length_x"), py::arg("length_y"),
             py::arg("length_z"), py::arg("viscosity"))
        .def("set_velocity", &StaggeredNavierStokesSolver3D::setVelocity, py::arg("i"), py::arg("j"),
             py::arg("k"), py::arg("u"), py::arg("v"), py::arg("w"))
        .def("stable_time_step", &StaggeredNavierStokesSolver3D::stableTimeStep, py::arg("velocity_scale"))
        .def("step", &StaggeredNavierStokesSolver3D::step, py::arg("dt"))
        .def("u", &StaggeredNavierStokesSolver3D::u)
        .def("v", &StaggeredNavierStokesSolver3D::v)
        .def("w", &StaggeredNavierStokesSolver3D::w)
        .def("pressure", &StaggeredNavierStokesSolver3D::pressure)
        .def("time", &StaggeredNavierStokesSolver3D::time)
        .def("max_divergence", &StaggeredNavierStokesSolver3D::maxDivergence);

    py::class_<StaggeredLidDrivenCavitySolver3D>(m, "StaggeredLidDrivenCavitySolver3D")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, double, double, double, double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("length_x"), py::arg("length_y"),
             py::arg("length_z"), py::arg("viscosity"), py::arg("lid_velocity"))
        .def("stable_time_step", &StaggeredLidDrivenCavitySolver3D::stableTimeStep)
        .def("step", &StaggeredLidDrivenCavitySolver3D::step, py::arg("dt"))
        .def("u", &StaggeredLidDrivenCavitySolver3D::u)
        .def("v", &StaggeredLidDrivenCavitySolver3D::v)
        .def("w", &StaggeredLidDrivenCavitySolver3D::w)
        .def("pressure", &StaggeredLidDrivenCavitySolver3D::pressure)
        .def("time", &StaggeredLidDrivenCavitySolver3D::time)
        .def("max_divergence", &StaggeredLidDrivenCavitySolver3D::maxDivergence)
        .def("load_state", &StaggeredLidDrivenCavitySolver3D::loadState, py::arg("u"), py::arg("v"), py::arg("w"),
             py::arg("p"), py::arg("time"));

    py::class_<MixingLengthLidDrivenCavitySolver3D>(m, "MixingLengthLidDrivenCavitySolver3D")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, double, double, double, double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("length_x"), py::arg("length_y"),
             py::arg("length_z"), py::arg("viscosity"), py::arg("lid_velocity"))
        .def("stable_time_step", &MixingLengthLidDrivenCavitySolver3D::stableTimeStep)
        .def("step", &MixingLengthLidDrivenCavitySolver3D::step, py::arg("dt"))
        .def("u", &MixingLengthLidDrivenCavitySolver3D::u)
        .def("v", &MixingLengthLidDrivenCavitySolver3D::v)
        .def("w", &MixingLengthLidDrivenCavitySolver3D::w)
        .def("pressure", &MixingLengthLidDrivenCavitySolver3D::pressure)
        .def("eddy_viscosity", &MixingLengthLidDrivenCavitySolver3D::eddyViscosity)
        .def("time", &MixingLengthLidDrivenCavitySolver3D::time)
        .def("max_divergence", &MixingLengthLidDrivenCavitySolver3D::maxDivergence);

    py::class_<SmagorinskyLesLidDrivenCavitySolver3D>(m, "SmagorinskyLesLidDrivenCavitySolver3D")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, double, double, double, double,
                      double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("length_x"), py::arg("length_y"),
             py::arg("length_z"), py::arg("viscosity"), py::arg("lid_velocity"),
             py::arg("smagorinsky_constant") = 0.17)
        .def("stable_time_step", &SmagorinskyLesLidDrivenCavitySolver3D::stableTimeStep)
        .def("step", &SmagorinskyLesLidDrivenCavitySolver3D::step, py::arg("dt"))
        .def("u", &SmagorinskyLesLidDrivenCavitySolver3D::u)
        .def("v", &SmagorinskyLesLidDrivenCavitySolver3D::v)
        .def("w", &SmagorinskyLesLidDrivenCavitySolver3D::w)
        .def("pressure", &SmagorinskyLesLidDrivenCavitySolver3D::pressure)
        .def("subgrid_viscosity", &SmagorinskyLesLidDrivenCavitySolver3D::subgridViscosity)
        .def("filter_width", &SmagorinskyLesLidDrivenCavitySolver3D::filterWidth)
        .def("time", &SmagorinskyLesLidDrivenCavitySolver3D::time)
        .def("max_divergence", &SmagorinskyLesLidDrivenCavitySolver3D::maxDivergence);

    py::class_<DesSstLidDrivenCavitySolver3D>(m, "DesSstLidDrivenCavitySolver3D")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, double, double, double, double, double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("length_x"), py::arg("length_y"),
             py::arg("length_z"), py::arg("viscosity"), py::arg("lid_velocity"), py::arg("c_des") = 0.61)
        .def("stable_time_step", &DesSstLidDrivenCavitySolver3D::stableTimeStep)
        .def("step", &DesSstLidDrivenCavitySolver3D::step, py::arg("dt"))
        .def("u", &DesSstLidDrivenCavitySolver3D::u)
        .def("v", &DesSstLidDrivenCavitySolver3D::v)
        .def("w", &DesSstLidDrivenCavitySolver3D::w)
        .def("pressure", &DesSstLidDrivenCavitySolver3D::pressure)
        .def("k", &DesSstLidDrivenCavitySolver3D::k)
        .def("omega", &DesSstLidDrivenCavitySolver3D::omega)
        .def("eddy_viscosity", &DesSstLidDrivenCavitySolver3D::eddyViscosity)
        .def("des_factor", &DesSstLidDrivenCavitySolver3D::desFactor)
        .def("filter_width", &DesSstLidDrivenCavitySolver3D::filterWidth)
        .def("time", &DesSstLidDrivenCavitySolver3D::time)
        .def("max_divergence", &DesSstLidDrivenCavitySolver3D::maxDivergence);

    py::class_<ImplicitConvectionDiffusionSolver1D> implicitConvectionDiffusion(
        m, "ImplicitConvectionDiffusionSolver1D");
    py::enum_<ImplicitConvectionDiffusionSolver1D::Preconditioner>(implicitConvectionDiffusion,
                                                                    "Preconditioner")
        .value("NONE", ImplicitConvectionDiffusionSolver1D::Preconditioner::None)
        .value("JACOBI", ImplicitConvectionDiffusionSolver1D::Preconditioner::Jacobi)
        .value("INCOMPLETE_LU", ImplicitConvectionDiffusionSolver1D::Preconditioner::IncompleteLU);
    implicitConvectionDiffusion
        .def(py::init<std::size_t, double, double, double, double, double, double>(), py::arg("nx"),
             py::arg("length"), py::arg("velocity"), py::arg("diffusivity"), py::arg("source"),
             py::arg("left_value"), py::arg("right_value"))
        .def("set_velocity_field", &ImplicitConvectionDiffusionSolver1D::setVelocityField,
             py::arg("velocity_per_cell"))
        .def("set_diffusivity_field", &ImplicitConvectionDiffusionSolver1D::setDiffusivityField,
             py::arg("diffusivity_per_cell"))
        .def("set_preconditioner", &ImplicitConvectionDiffusionSolver1D::setPreconditioner,
             py::arg("preconditioner"))
        .def("solve_gauss_seidel", &ImplicitConvectionDiffusionSolver1D::solveGaussSeidel,
             py::arg("max_iterations") = 20000, py::arg("tolerance") = 1e-12)
        .def("solve_bicgstab", &ImplicitConvectionDiffusionSolver1D::solveBiCGStab,
             py::arg("max_iterations") = 500, py::arg("tolerance") = 1e-10)
        .def("solve_gmres", &ImplicitConvectionDiffusionSolver1D::solveGmres, py::arg("restart") = 30,
             py::arg("max_iterations") = 500, py::arg("tolerance") = 1e-10)
        .def("residual_history", &ImplicitConvectionDiffusionSolver1D::residualHistory)
        .def("value", &ImplicitConvectionDiffusionSolver1D::value)
        .def("cell_center_x", &ImplicitConvectionDiffusionSolver1D::cellCenterX);

    py::enum_<MultigridPoissonSolver2D::Face>(m, "MultigridFace")
        .value("X_MIN", MultigridPoissonSolver2D::Face::XMin)
        .value("X_MAX", MultigridPoissonSolver2D::Face::XMax)
        .value("Y_MIN", MultigridPoissonSolver2D::Face::YMin)
        .value("Y_MAX", MultigridPoissonSolver2D::Face::YMax);

    py::enum_<ImplicitConvectionDiffusionSolver1D::ConvectionScheme>(implicitConvectionDiffusion,
                                                                     "ConvectionScheme")
        .value("FIRST_ORDER_UPWIND",
               ImplicitConvectionDiffusionSolver1D::ConvectionScheme::FirstOrderUpwind)
        .value("CENTRAL", ImplicitConvectionDiffusionSolver1D::ConvectionScheme::Central)
        .value("LIMITED_LINEAR_UPWIND",
               ImplicitConvectionDiffusionSolver1D::ConvectionScheme::LimitedLinearUpwind);
    implicitConvectionDiffusion
        .def("set_convection_scheme", &ImplicitConvectionDiffusionSolver1D::setConvectionScheme,
             py::arg("scheme"))
        // Zero is a theorem here, not a tolerance: the equation obeys a
        // maximum principle, so any overshoot is the scheme failing.
        .def("max_boundedness_violation",
             &ImplicitConvectionDiffusionSolver1D::maxBoundednessViolation)
        .def("max_cell_peclet", &ImplicitConvectionDiffusionSolver1D::maxCellPeclet);

    py::class_<MultigridPoissonSolver2D>(m, "MultigridPoissonSolver2D")
        .def(py::init<std::size_t, std::size_t, double, double>(), py::arg("nx"), py::arg("ny"),
             py::arg("length_x"), py::arg("length_y"))
        .def("set_boundary_value", &MultigridPoissonSolver2D::setBoundaryValue, py::arg("face"),
             py::arg("value"))
        .def("set_source_term", &MultigridPoissonSolver2D::setSourceTerm, py::arg("source"))
        .def("solve", &MultigridPoissonSolver2D::solve, py::arg("max_v_cycles") = 100,
             py::arg("tolerance") = 1e-9, py::arg("pre_sweeps") = 2, py::arg("post_sweeps") = 2)
        .def("value", &MultigridPoissonSolver2D::value);

    py::class_<KOmegaSSTLidDrivenCavitySolver2D>(m, "KOmegaSSTLidDrivenCavitySolver2D")
        .def(py::init<std::size_t, std::size_t, double, double, double, double>(), py::arg("nx"),
             py::arg("ny"), py::arg("length_x"), py::arg("length_y"), py::arg("viscosity"),
             py::arg("lid_velocity"))
        .def("stable_time_step", &KOmegaSSTLidDrivenCavitySolver2D::stableTimeStep)
        .def("step", &KOmegaSSTLidDrivenCavitySolver2D::step, py::arg("dt"))
        .def("u", &KOmegaSSTLidDrivenCavitySolver2D::u)
        .def("v", &KOmegaSSTLidDrivenCavitySolver2D::v)
        .def("pressure", &KOmegaSSTLidDrivenCavitySolver2D::pressure)
        .def("k", &KOmegaSSTLidDrivenCavitySolver2D::k)
        .def("omega", &KOmegaSSTLidDrivenCavitySolver2D::omega)
        .def("eddy_viscosity", &KOmegaSSTLidDrivenCavitySolver2D::eddyViscosity)
        .def("time", &KOmegaSSTLidDrivenCavitySolver2D::time)
        .def("max_divergence", &KOmegaSSTLidDrivenCavitySolver2D::maxDivergence);

    py::class_<KOmegaSSTLidDrivenCavitySolver3D>(m, "KOmegaSSTLidDrivenCavitySolver3D")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, double, double, double, double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("length_x"), py::arg("length_y"),
             py::arg("length_z"), py::arg("viscosity"), py::arg("lid_velocity"))
        .def("stable_time_step", &KOmegaSSTLidDrivenCavitySolver3D::stableTimeStep)
        .def("step", &KOmegaSSTLidDrivenCavitySolver3D::step, py::arg("dt"))
        .def("u", &KOmegaSSTLidDrivenCavitySolver3D::u)
        .def("v", &KOmegaSSTLidDrivenCavitySolver3D::v)
        .def("w", &KOmegaSSTLidDrivenCavitySolver3D::w)
        .def("pressure", &KOmegaSSTLidDrivenCavitySolver3D::pressure)
        .def("k", &KOmegaSSTLidDrivenCavitySolver3D::k)
        .def("omega", &KOmegaSSTLidDrivenCavitySolver3D::omega)
        .def("eddy_viscosity", &KOmegaSSTLidDrivenCavitySolver3D::eddyViscosity)
        .def("time", &KOmegaSSTLidDrivenCavitySolver3D::time)
        .def("max_divergence", &KOmegaSSTLidDrivenCavitySolver3D::maxDivergence);

    py::class_<KEpsilonLidDrivenCavitySolver2D>(m, "KEpsilonLidDrivenCavitySolver2D")
        .def(py::init<std::size_t, std::size_t, double, double, double, double>(), py::arg("nx"),
             py::arg("ny"), py::arg("length_x"), py::arg("length_y"), py::arg("viscosity"),
             py::arg("lid_velocity"))
        .def("stable_time_step", &KEpsilonLidDrivenCavitySolver2D::stableTimeStep)
        .def("step", &KEpsilonLidDrivenCavitySolver2D::step, py::arg("dt"))
        .def("u", &KEpsilonLidDrivenCavitySolver2D::u)
        .def("v", &KEpsilonLidDrivenCavitySolver2D::v)
        .def("pressure", &KEpsilonLidDrivenCavitySolver2D::pressure)
        .def("k", &KEpsilonLidDrivenCavitySolver2D::k)
        .def("epsilon", &KEpsilonLidDrivenCavitySolver2D::epsilon)
        .def("eddy_viscosity", &KEpsilonLidDrivenCavitySolver2D::eddyViscosity)
        .def("time", &KEpsilonLidDrivenCavitySolver2D::time)
        .def("max_divergence", &KEpsilonLidDrivenCavitySolver2D::maxDivergence);

    py::class_<KEpsilonLidDrivenCavitySolver3D>(m, "KEpsilonLidDrivenCavitySolver3D")
        .def(py::init<std::size_t, std::size_t, std::size_t, double, double, double, double, double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("length_x"), py::arg("length_y"),
             py::arg("length_z"), py::arg("viscosity"), py::arg("lid_velocity"))
        .def("stable_time_step", &KEpsilonLidDrivenCavitySolver3D::stableTimeStep)
        .def("step", &KEpsilonLidDrivenCavitySolver3D::step, py::arg("dt"))
        .def("u", &KEpsilonLidDrivenCavitySolver3D::u)
        .def("v", &KEpsilonLidDrivenCavitySolver3D::v)
        .def("w", &KEpsilonLidDrivenCavitySolver3D::w)
        .def("pressure", &KEpsilonLidDrivenCavitySolver3D::pressure)
        .def("k", &KEpsilonLidDrivenCavitySolver3D::k)
        .def("epsilon", &KEpsilonLidDrivenCavitySolver3D::epsilon)
        .def("eddy_viscosity", &KEpsilonLidDrivenCavitySolver3D::eddyViscosity)
        .def("time", &KEpsilonLidDrivenCavitySolver3D::time)
        .def("max_divergence", &KEpsilonLidDrivenCavitySolver3D::maxDivergence);

    py::class_<MixingLengthLidDrivenCavitySolver2D>(m, "MixingLengthLidDrivenCavitySolver2D")
        .def(py::init<std::size_t, std::size_t, double, double, double, double>(), py::arg("nx"),
             py::arg("ny"), py::arg("length_x"), py::arg("length_y"), py::arg("viscosity"),
             py::arg("lid_velocity"))
        .def("stable_time_step", &MixingLengthLidDrivenCavitySolver2D::stableTimeStep)
        .def("step", &MixingLengthLidDrivenCavitySolver2D::step, py::arg("dt"))
        .def("u", &MixingLengthLidDrivenCavitySolver2D::u)
        .def("v", &MixingLengthLidDrivenCavitySolver2D::v)
        .def("pressure", &MixingLengthLidDrivenCavitySolver2D::pressure)
        .def("eddy_viscosity", &MixingLengthLidDrivenCavitySolver2D::eddyViscosity)
        .def("time", &MixingLengthLidDrivenCavitySolver2D::time)
        .def("max_divergence", &MixingLengthLidDrivenCavitySolver2D::maxDivergence);

    py::class_<MixingLengthChannelFlowSolver1D>(m, "MixingLengthChannelFlowSolver1D")
        .def(py::init<std::size_t, double, double, double>(), py::arg("ny"), py::arg("height"),
             py::arg("kinematic_viscosity"), py::arg("source"))
        .def("solve", &MixingLengthChannelFlowSolver1D::solve, py::arg("max_outer_iterations") = 500,
             py::arg("tolerance") = 1e-10)
        .def("u", &MixingLengthChannelFlowSolver1D::u)
        .def("eddy_viscosity", &MixingLengthChannelFlowSolver1D::eddyViscosity)
        .def("wall_distance", &MixingLengthChannelFlowSolver1D::wallDistance)
        .def("cell_center_y", &MixingLengthChannelFlowSolver1D::cellCenterY)
        .def("friction_velocity", &MixingLengthChannelFlowSolver1D::frictionVelocity);

    py::class_<KEpsilonChannelFlowSolver1D>(m, "KEpsilonChannelFlowSolver1D")
        .def(py::init<std::size_t, double, double, double>(), py::arg("ny"), py::arg("height"),
             py::arg("kinematic_viscosity"), py::arg("source"))
        .def("solve", &KEpsilonChannelFlowSolver1D::solve, py::arg("max_outer_iterations") = 2000,
             py::arg("tolerance") = 1e-10)
        .def("u", &KEpsilonChannelFlowSolver1D::u)
        .def("k", &KEpsilonChannelFlowSolver1D::k)
        .def("epsilon", &KEpsilonChannelFlowSolver1D::epsilon)
        .def("eddy_viscosity", &KEpsilonChannelFlowSolver1D::eddyViscosity)
        .def("wall_distance", &KEpsilonChannelFlowSolver1D::wallDistance)
        .def("cell_center_y", &KEpsilonChannelFlowSolver1D::cellCenterY)
        .def("friction_velocity", &KEpsilonChannelFlowSolver1D::frictionVelocity);

    py::class_<KOmegaSSTChannelFlowSolver1D>(m, "KOmegaSSTChannelFlowSolver1D")
        .def(py::init<std::size_t, double, double, double>(), py::arg("ny"), py::arg("height"),
             py::arg("kinematic_viscosity"), py::arg("source"))
        .def("solve", &KOmegaSSTChannelFlowSolver1D::solve, py::arg("max_outer_iterations") = 2000,
             py::arg("tolerance") = 1e-10)
        .def("u", &KOmegaSSTChannelFlowSolver1D::u)
        .def("k", &KOmegaSSTChannelFlowSolver1D::k)
        .def("omega", &KOmegaSSTChannelFlowSolver1D::omega)
        .def("eddy_viscosity", &KOmegaSSTChannelFlowSolver1D::eddyViscosity)
        .def("wall_distance", &KOmegaSSTChannelFlowSolver1D::wallDistance)
        .def("cell_center_y", &KOmegaSSTChannelFlowSolver1D::cellCenterY)
        .def("friction_velocity", &KOmegaSSTChannelFlowSolver1D::frictionVelocity);

    // -- Unstructured FVM on tetrahedra (ROADMAP Fase 2 and 3) -------------
    //
    // These were the one layer of the engine without bindings, which is not
    // a cosmetic gap: every experiment with them meant writing and compiling
    // C++, which is exactly the friction the Python layer exists to remove.
    // The investigation that found the 13.2% mass imbalance would have been
    // minutes here instead of several build cycles -- that is the argument,
    // recorded as DIVIDA_TECNICA.md 5.2.
    //
    // Both classes take the mesh by const reference and keep a pointer to
    // it, hence keep_alive: letting Python collect the mesh while a solver
    // still points at it is a use-after-free, and the mesh is the natural
    // thing for a script to drop after constructing the solver.
    //
    // cell_count / max_non_orthogonality / deficient_stencil_count are
    // inherited from UnstructuredFvmBase and re-declared on each class
    // rather than binding the base itself, which is not constructible or
    // destructible from outside.

    py::class_<UnstructuredDiffusionSolver>(m, "UnstructuredDiffusionSolver")
        .def(py::init<const aether::mesh::TetrahedralMesh&>(), py::arg("mesh"), py::keep_alive<1, 2>())
        // `selector` is called once per boundary face with that face's
        // centroid: "the face on the x = 0 plane" without needing to know
        // how the mesh generator numbered anything. Every face starts
        // insulated; a later call overrides an earlier one.
        .def("set_dirichlet_boundary",
             py::overload_cast<const std::function<bool(const aether::core::Vector3&)>&, double>(
                 &UnstructuredDiffusionSolver::setDirichletBoundary),
             py::arg("selector"), py::arg("value"))
        // Value varying over the boundary, for data that is not piecewise
        // constant -- a manufactured solution above all. Registered second so
        // a plain number still binds to the constant overload.
        .def("set_dirichlet_boundary",
             py::overload_cast<const std::function<bool(const aether::core::Vector3&)>&,
                                const std::function<double(const aether::core::Vector3&)>&>(
                 &UnstructuredDiffusionSolver::setDirichletBoundary),
             py::arg("selector"), py::arg("value"))
        // Volumetric source: turns Laplace into Poisson, which is what the
        // method of manufactured solutions needs -- pick any smooth phi, set
        // S = -laplacian(phi), and the discretization error is all that is
        // left to measure.
        .def("set_source_term", &UnstructuredDiffusionSolver::setSourceTerm, py::arg("source"))
        .def("solve_conjugate_gradient", &UnstructuredDiffusionSolver::solveConjugateGradient,
             py::arg("max_iterations") = 20000, py::arg("tolerance") = 1e-10,
             py::arg("max_outer_sweeps") = 50)
        .def("value", &UnstructuredDiffusionSolver::value, py::arg("cell"))
        .def("cell_gradients", &UnstructuredDiffusionSolver::cellGradients)
        // Largest per-cell change across the final outer sweep: what
        // separates a solve that settled from one that hit the sweep cap.
        .def("last_outer_change", &UnstructuredDiffusionSolver::lastOuterChange)
        .def("cell_count", &UnstructuredDiffusionSolver::cellCount)
        .def("max_non_orthogonality", &UnstructuredDiffusionSolver::maxNonOrthogonality)
        .def("deficient_stencil_count", &UnstructuredDiffusionSolver::deficientStencilCount);

    // Steady convection-diffusion of a passive scalar carried by a
    // *prescribed* velocity. Exists to measure the convection scheme against
    // a known answer -- the Navier-Stokes solver cannot do that, because
    // there the velocity is what is being solved for.
    py::class_<UnstructuredScalarTransportSolver> transport(m, "UnstructuredScalarTransportSolver");
    py::enum_<UnstructuredScalarTransportSolver::ConvectionScheme>(transport, "ConvectionScheme")
        .value("FIRST_ORDER_UPWIND",
               UnstructuredScalarTransportSolver::ConvectionScheme::FirstOrderUpwind)
        .value("LIMITED_LINEAR_UPWIND",
               UnstructuredScalarTransportSolver::ConvectionScheme::LimitedLinearUpwind);
    transport
        .def(py::init<const aether::mesh::TetrahedralMesh&, double,
                      std::function<aether::core::Vector3(const aether::core::Vector3&)>,
                      UnstructuredScalarTransportSolver::ConvectionScheme>(),
             py::arg("mesh"), py::arg("diffusivity"), py::arg("velocity"),
             py::arg("scheme") =
                 UnstructuredScalarTransportSolver::ConvectionScheme::LimitedLinearUpwind,
             py::keep_alive<1, 2>())
        .def("set_dirichlet_boundary", &UnstructuredScalarTransportSolver::setDirichletBoundary,
             py::arg("selector"), py::arg("value"))
        .def("set_source_term", &UnstructuredScalarTransportSolver::setSourceTerm, py::arg("source"))
        .def("solve_steady", &UnstructuredScalarTransportSolver::solveSteady,
             py::arg("tolerance") = 1e-12, py::arg("max_steps") = 500000)
        .def("stable_time_step", &UnstructuredScalarTransportSolver::stableTimeStep)
        .def("value", &UnstructuredScalarTransportSolver::value, py::arg("cell"))
        .def("last_change", &UnstructuredScalarTransportSolver::lastChange)
        // How far into the convection-dominated regime this case sits, which
        // is the only regime where the error measures the convection scheme.
        .def("max_cell_peclet", &UnstructuredScalarTransportSolver::maxCellPeclet)
        .def("cell_count", &UnstructuredScalarTransportSolver::cellCount)
        .def("max_non_orthogonality", &UnstructuredScalarTransportSolver::maxNonOrthogonality)
        .def("deficient_stencil_count", &UnstructuredScalarTransportSolver::deficientStencilCount);

    py::class_<UnstructuredCavitySolver3D> unstructuredCavity(m, "UnstructuredCavitySolver3D");
    py::enum_<UnstructuredCavitySolver3D::EnergyModel>(unstructuredCavity, "EnergyModel")
        .value("NONE", UnstructuredCavitySolver3D::EnergyModel::None)
        .value("PASSIVE", UnstructuredCavitySolver3D::EnergyModel::Passive)
        .value("BOUSSINESQ", UnstructuredCavitySolver3D::EnergyModel::Boussinesq);
    py::enum_<UnstructuredCavitySolver3D::TurbulenceModel>(unstructuredCavity, "TurbulenceModel")
        .value("NONE", UnstructuredCavitySolver3D::TurbulenceModel::None)
        .value("MIXING_LENGTH", UnstructuredCavitySolver3D::TurbulenceModel::MixingLength);
    unstructuredCavity
        // `wall_velocity(position)` gives the prescribed velocity at a
        // boundary face centroid -- zero on a solid wall, the lid's speed on
        // a moving one, the inlet profile on an inlet (an inlet needs
        // nothing new: it is a wall with a non-zero prescribed velocity).
        // `is_outlet(position)` marks faces where fluid may leave; pass None
        // for a closed domain, where every boundary is a solid wall.
        // `pressure_correctors` is exposed because how many the projection
        // needs is a property of the *mesh*: four leaves 1e-04 of the inflow
        // unaccounted on a jittered mesh that sixty-four closes to 1e-14.
        // Pair it with last_pressure_change() rather than guessing.
        //
        // `turbulence` selects the closure; the convection scheme is left at
        // the class default rather than exposed alongside it, because the
        // measurement that settled that choice (DIVIDA_TECNICA.md 3.1) is
        // not one a caller should be invited to re-litigate per call.
        .def(py::init([](const aether::mesh::TetrahedralMesh& mesh, double viscosity,
                          std::function<aether::core::Vector3(const aether::core::Vector3&)> wallVelocity,
                          std::function<bool(const aether::core::Vector3&)> isOutlet,
                          double outletPressure, std::size_t pressureCorrectors,
                          UnstructuredCavitySolver3D::TurbulenceModel turbulence) {
                 return std::make_unique<UnstructuredCavitySolver3D>(
                     mesh, viscosity, std::move(wallVelocity), std::move(isOutlet), outletPressure,
                     pressureCorrectors,
                     UnstructuredCavitySolver3D::ConvectionScheme::LimitedLinearUpwind, turbulence);
             }),
             py::arg("mesh"), py::arg("viscosity"), py::arg("wall_velocity"),
             py::arg("is_outlet") = py::none(), py::arg("outlet_pressure") = 0.0,
             py::arg("pressure_correctors") = UnstructuredCavitySolver3D::kDefaultPressureCorrectors,
             py::arg("turbulence") = UnstructuredCavitySolver3D::TurbulenceModel::None,
             py::keep_alive<1, 2>())
        .def("step", &UnstructuredCavitySolver3D::step, py::arg("dt"))
        // Measurement instrument: one step with pieces switched off, so the
        // step operator can be taken apart. See the header -- a step without
        // the projection does not conserve mass, and is not meant to.
        .def("step_with", [](UnstructuredCavitySolver3D& self, double dt, bool convection,
                              bool viscous, bool projection) {
                 UnstructuredCavitySolver3D::StepParts parts;
                 parts.convection = convection;
                 parts.viscous = viscous;
                 parts.projection = projection;
                 self.stepWith(dt, parts);
             },
             py::arg("dt"), py::arg("convection") = true, py::arg("viscous") = true,
             py::arg("projection") = true)
        // Imposes an initial field. Checkpointing is one use; the other is
        // measuring the step operator itself, which is what
        // DIVIDA_TECNICA.md 4.3 was blocked on.
        .def("load_state", &UnstructuredCavitySolver3D::loadState, py::arg("velocity"),
             py::arg("pressure"), py::arg("time") = 0.0)
        .def("stable_time_step", &UnstructuredCavitySolver3D::stableTimeStep)
        .def("velocity", &UnstructuredCavitySolver3D::velocity, py::arg("cell"))
        .def("pressure", &UnstructuredCavitySolver3D::pressure, py::arg("cell"))
        .def("time", &UnstructuredCavitySolver3D::time)
        .def("cell_count", &UnstructuredCavitySolver3D::cellCount)
        // The quantity the projection actually drives to zero, measured on
        // faces rather than as a cell-centred difference -- see the class
        // comment for why that distinction has bitten this project twice.
        .def("max_face_divergence", &UnstructuredCavitySolver3D::maxFaceDivergence)
        // Net mass flux across every non-wall boundary face. A per-cell
        // divergence can be zero everywhere while the domain as a whole
        // gains or loses mass, which is what this catches.
        .def("net_boundary_flux", &UnstructuredCavitySolver3D::netBoundaryFlux)
        .def("max_non_orthogonality", &UnstructuredCavitySolver3D::maxNonOrthogonality)
        // Whether the projection converged or merely ran out of correctors.
        .def("last_pressure_change", &UnstructuredCavitySolver3D::lastPressureChange)
        // Below 1.0 means this mesh forced the correction to be damped.
        .def("pressure_relaxation", &UnstructuredCavitySolver3D::pressureRelaxation)
        // The GMRES coupling correction from DIVIDA_TECNICA.md 4.3: how many
        // matrix-vector products the last step's correction used (0 means
        // the mesh did not need one), and the residual before/after.
        .def("last_coupling_iterations", &UnstructuredCavitySolver3D::lastCouplingIterations)
        .def("last_coupling_residual_before", &UnstructuredCavitySolver3D::lastCouplingResidualBefore)
        .def("last_coupling_residual_after", &UnstructuredCavitySolver3D::lastCouplingResidualAfter)
        // How many cells fall back to Green-Gauss because their
        // least-squares stencil is rank-deficient. Not a curiosity: it was
        // 7% of the cavity and 12% of the channel, and those cells used to
        // receive a silent zero gradient instead.
        .def("deficient_stencil_count", &UnstructuredCavitySolver3D::deficientStencilCount)
        // Turbulent viscosity per cell (identically zero without a closure,
        // and zero at any solid wall even with one), and the wall distance
        // the mixing length is built from.
        .def("eddy_viscosity", &UnstructuredCavitySolver3D::eddyViscosity, py::arg("cell"))
        .def("wall_distance", &UnstructuredCavitySolver3D::wallDistance, py::arg("cell"))
        // Temperature transport (ROADMAP Fase 6.1). PASSIVE is forced
        // convection -- the flow carries the heat and the heat does not move
        // the flow; BOUSSINESQ adds the buoyancy feedback, i.e. natural
        // convection. A boundary face no set_wall_temperature() call ever
        // selects stays adiabatic.
        .def("enable_energy", &UnstructuredCavitySolver3D::enableEnergy, py::arg("model"),
             py::arg("thermal_diffusivity"), py::arg("reference_temperature") = 0.0,
             py::arg("thermal_expansion") = 0.0,
             py::arg("gravity") = aether::core::Vector3(0.0, 0.0, -9.81))
        .def("set_wall_temperature", &UnstructuredCavitySolver3D::setWallTemperature,
             py::arg("selector"), py::arg("value"))
        .def("temperature", &UnstructuredCavitySolver3D::temperature, py::arg("cell"))
        .def("set_temperature", &UnstructuredCavitySolver3D::setTemperature, py::arg("cell"),
             py::arg("value"));
}
