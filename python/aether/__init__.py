"""High-level Python orchestration layer for the Aether CFD engine.

The numerical core, geometry, mesh and solver layers live in C++ under
engine/{core,geometry,mesh,solver}, exposed here through the compiled
aether_{core,geometry,mesh,solver}_py extensions built from bindings/python.

Import order matters: aether_core_py must load first because it registers
aether::core::Vector3 with pybind11, and the geometry/mesh/solver bindings
return or accept that type (and StructuredGrid3D, for the solver) without
re-registering it themselves.

engine/gpu (Module 10, CUDA) is the one layer imported *optionally*: the
aether_gpu_py extension only exists at all on a machine where CMake found
the CUDA Toolkit (see bindings/python/CMakeLists.txt's own TARGET guard),
so a missing extension here is expected, not an error -- check the
GPU_AVAILABLE flag below rather than assuming PoissonOperatorCuda works.
"""



def _ensure_extensions_importable() -> None:
    """Puts the compiled extensions on sys.path when they are not already
    there, by looking where this repository's own build actually writes
    them.

    **Why this exists.** Until this was added, every entry point into the
    package -- a script, a test, an interactive session -- had to set
    PYTHONPATH to the build output directory by hand, and that directory
    differs by generator: a multi-config generator (Visual Studio) adds a
    Release/ or Debug/ level that a single-config one does not. It is the
    first thing anyone hits and the least interesting, so the package
    resolves it rather than documenting it.

    An extension that is *already* importable is left alone and this
    returns immediately: an explicit PYTHONPATH, a virtualenv install or a
    future wheel all mean the caller has already said where the extensions
    are, and overriding that would be worse than the friction it removes.
    """
    import importlib.util
    import pathlib
    import sys

    if importlib.util.find_spec("aether_core_py") is not None:
        return

    # python/aether/__init__.py -> python/aether -> python -> repository root
    repo_root = pathlib.Path(__file__).resolve().parent.parent.parent
    build = repo_root / "build" / "bindings" / "python"
    for candidate in (build / "Release", build / "Debug", build):
        if candidate.is_dir() and any(candidate.glob("aether_core_py*")):
            sys.path.insert(0, str(candidate))
            return


_ensure_extensions_importable()

try:
    from aether_core_py import Mesh, ScalarField, Tensor3x3, Vector3, VectorField
    from aether_geometry_py import TriangleMesh, load_obj, load_stl, save_obj, save_stl_binary
    from aether_mesh_py import (
        DelaunayTetrahedralization3D,
        DelaunayTriangulation2D,
        PolygonTriangulation2D,
        StructuredGrid3D,
        TetrahedralFace,
        TetrahedralMesh,
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
        UnstructuredCavitySolver3D,
        UnstructuredDiffusionSolver,
        UnstructuredScalarTransportSolver,
    )
    from aether_postprocessing_py import (
        Segment2D,
        Streamline2D,
        Triangle3D,
        marching_cubes_3d,
        marching_squares_2d,
        write_tetrahedral_mesh_vtk,
    )
    from aether_persistence_py import (
        FieldArchive,
        HistoryEntry,
        ProjectHistory,
        load_grid,
        load_tetrahedral_mesh,
        save_grid,
        save_tetrahedral_mesh,
    )
    from aether_optimization_py import NelderMead, OptimizationResult
    from aether_analysis_py import FieldStatistics, checkerboard_index, compute_statistics, max_courant_number, summarize_field
    from aether_ml_py import MultiLayerPerceptron
    from aether_plugin_py import DiagnosticInfo, PluginHost
    from aether.pipeline import (
        ConservationReport,
        FlowDomain,
        MeshGenerationError,
        RunReport,
        StabilityReport,
        check_closed_domain_conservation,
        classify_boundary_face,
        driving_wall_velocity,
        export_result_vtk,
        measure_mesh_stability,
        mesh_flow_around_object,
        run_to_steady_state,
    )
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "aether_*_py extensions not found. Build the C++ core first:\n"
        "  python build.py\n"
        "(or by hand: cmake -S . -B build && cmake --build build --config Release).\n"
        "The build output is found automatically when it sits in this repository's\n"
        "own build/ directory; set PYTHONPATH to it explicitly if you built elsewhere."
    ) from exc

# Separate, non-fatal try/except: unlike every extension above, not
# building this one is a normal outcome (no CUDA Toolkit on this machine),
# not a broken build -- so its absence must not take down the rest of the
# package the way a missing core/geometry/mesh/solver/postprocessing
# extension does above.
try:
    from aether_gpu_py import PoissonOperatorCuda

    GPU_AVAILABLE = True
except ImportError:
    PoissonOperatorCuda = None
    GPU_AVAILABLE = False

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
    "TetrahedralMesh",
    "TetrahedralFace",
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
    "UnstructuredDiffusionSolver",
    "UnstructuredCavitySolver3D",
    "UnstructuredScalarTransportSolver",
    "Streamline2D",
    "Segment2D",
    "marching_squares_2d",
    "Triangle3D",
    "marching_cubes_3d",
    "FieldArchive",
    "HistoryEntry",
    "ProjectHistory",
    "save_grid",
    "load_grid",
    "save_tetrahedral_mesh",
    "load_tetrahedral_mesh",
    "NelderMead",
    "OptimizationResult",
    "FieldStatistics",
    "compute_statistics",
    "max_courant_number",
    "checkerboard_index",
    "summarize_field",
    "MultiLayerPerceptron",
    "PluginHost",
    "DiagnosticInfo",
    "PoissonOperatorCuda",
    "GPU_AVAILABLE",
    "mesh_flow_around_object",
    "classify_boundary_face",
    "driving_wall_velocity",
    "FlowDomain",
    "MeshGenerationError",
    "check_closed_domain_conservation",
    "ConservationReport",
    "measure_mesh_stability",
    "StabilityReport",
    "run_to_steady_state",
    "RunReport",
    "export_result_vtk",
    "write_tetrahedral_mesh_vtk",
]
