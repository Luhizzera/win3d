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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

    py::class_<LidDrivenCavitySolver2D>(m, "LidDrivenCavitySolver2D")
        .def(py::init<std::size_t, std::size_t, double, double, double, double>(), py::arg("nx"),
             py::arg("ny"), py::arg("length_x"), py::arg("length_y"), py::arg("viscosity"),
             py::arg("lid_velocity"))
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
}
