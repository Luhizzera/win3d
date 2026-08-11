"""High-level Python orchestration layer for the Aether CFD engine.

The numerical core, geometry, mesh and solver layers live in C++ under
engine/{core,geometry,mesh,solver}, exposed here through the compiled
aether_{core,geometry,mesh,solver}_py extensions built from bindings/python.

Import order matters: aether_core_py must load first because it registers
aether::core::Vector3 with pybind11, and the geometry/mesh/solver bindings
return or accept that type (and StructuredGrid3D, for the solver) without
re-registering it themselves.
"""

try:
    from aether_core_py import Mesh, ScalarField, Tensor3x3, Vector3, VectorField
    from aether_geometry_py import TriangleMesh, load_obj, load_stl, save_obj, save_stl_binary
    from aether_mesh_py import (
        DelaunayTetrahedralization3D,
        DelaunayTriangulation2D,
        PolygonTriangulation2D,
        StructuredGrid3D,
    )
    from aether_solver_py import (
        DesSstLidDrivenCavitySolver3D,
        DiffusionProblem,
        Face,
        KEpsilonChannelFlowSolver1D,
        KEpsilonLidDrivenCavitySolver2D,
        KEpsilonLidDrivenCavitySolver3D,
        KOmegaSSTChannelFlowSolver1D,
        KOmegaSSTLidDrivenCavitySolver2D,
        KOmegaSSTLidDrivenCavitySolver3D,
        ImplicitConvectionDiffusionSolver1D,
        LidDrivenCavitySolver2D,
        MixingLengthChannelFlowSolver1D,
        MixingLengthLidDrivenCavitySolver2D,
        MixingLengthLidDrivenCavitySolver3D,
        MultigridFace,
        MultigridPoissonSolver2D,
        SmagorinskyLesLidDrivenCavitySolver3D,
        StaggeredLidDrivenCavitySolver3D,
        StaggeredNavierStokesSolver3D,
        SteadyDiffusionSolver,
        TaylorGreenVortexSolver2D,
        TransientDiffusionSolver,
    )
    from aether_postprocessing_py import (
        Segment2D,
        Streamline2D,
        Triangle3D,
        marching_cubes_3d,
        marching_squares_2d,
    )
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "aether_*_py extensions not found. Build the C++ core first:\n"
        "  cmake -S . -B build\n"
        "  cmake --build build\n"
        "then make sure the build output directory is on PYTHONPATH."
    ) from exc

__all__ = [
    "Vector3",
    "Tensor3x3",
    "Mesh",
    "ScalarField",
    "VectorField",
    "TriangleMesh",
    "load_stl",
    "save_stl_binary",
    "load_obj",
    "save_obj",
    "StructuredGrid3D",
    "DelaunayTriangulation2D",
    "PolygonTriangulation2D",
    "DelaunayTetrahedralization3D",
    "DiffusionProblem",
    "Face",
    "SteadyDiffusionSolver",
    "TransientDiffusionSolver",
    "TaylorGreenVortexSolver2D",
    "LidDrivenCavitySolver2D",
    "StaggeredNavierStokesSolver3D",
    "MixingLengthChannelFlowSolver1D",
    "KEpsilonChannelFlowSolver1D",
    "KEpsilonLidDrivenCavitySolver2D",
    "KEpsilonLidDrivenCavitySolver3D",
    "KOmegaSSTLidDrivenCavitySolver2D",
    "KOmegaSSTLidDrivenCavitySolver3D",
    "KOmegaSSTChannelFlowSolver1D",
    "MixingLengthLidDrivenCavitySolver2D",
    "MixingLengthLidDrivenCavitySolver3D",
    "MultigridPoissonSolver2D",
    "MultigridFace",
    "StaggeredLidDrivenCavitySolver3D",
    "SmagorinskyLesLidDrivenCavitySolver3D",
    "DesSstLidDrivenCavitySolver3D",
    "ImplicitConvectionDiffusionSolver1D",
    "Streamline2D",
    "Segment2D",
    "marching_squares_2d",
    "Triangle3D",
    "marching_cubes_3d",
]
