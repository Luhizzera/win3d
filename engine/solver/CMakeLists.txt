add_library(aether_solver STATIC
    src/DiffusionProblem.cpp
    src/ExplicitTimeStep.cpp
    src/SteadyDiffusionSolver.cpp
    src/TransientDiffusionSolver.cpp
    src/TaylorGreenVortexSolver2D.cpp
    src/LidDrivenCavitySolver2D.cpp
    src/UnstructuredFvmBase.cpp
    src/UnstructuredDiffusionSolver.cpp
    src/UnstructuredScalarTransportSolver.cpp
    src/UnstructuredCavitySolver3D.cpp
    src/MixingLengthLidDrivenCavitySolver2D.cpp
    src/MultigridPoissonSolver2D.cpp
    src/StaggeredCavityBase3D.cpp
    src/StaggeredLidDrivenCavitySolver3D.cpp
    src/MixingLengthLidDrivenCavitySolver3D.cpp
    src/KEpsilonLidDrivenCavitySolver3D.cpp
    src/KOmegaSSTLidDrivenCavitySolver3D.cpp
    src/SmagorinskyLesLidDrivenCavitySolver3D.cpp
    src/DesSstLidDrivenCavitySolver3D.cpp
    src/ImplicitConvectionDiffusionSolver1D.cpp
    src/StaggeredNavierStokesSolver3D.cpp
    src/MixingLengthChannelFlowSolver1D.cpp
    src/KEpsilonChannelFlowSolver1D.cpp
    src/KEpsilonLidDrivenCavitySolver2D.cpp
    src/KOmegaSSTChannelFlowSolver1D.cpp
    src/KOmegaSSTLidDrivenCavitySolver2D.cpp
    src/CompressibleEulerSolver1D.cpp
    src/LevelSetAdvectionSolver2D.cpp
)

target_include_directories(aether_solver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(aether_solver PUBLIC aether_mesh)

# Module 10 (GPU), opt-in and detected rather than required -- same
# if(TARGET aether_gpu) conditional bindings/python/CMakeLists.txt already
# uses for aether_gpu_py, correct regardless of AETHER_ENABLE_CUDA/
# CMAKE_CUDA_COMPILER scoping details between the two CMakeLists.txt files.
# AETHER_HAVE_GPU is PRIVATE: only StaggeredCavityBase3D.cpp's own #ifdef
# ever sees it (the public header never references it), so it cannot cause
# an ABI mismatch between translation units that include the header.
if(TARGET aether_gpu)
    target_link_libraries(aether_solver PRIVATE aether_gpu)
    target_compile_definitions(aether_solver PRIVATE AETHER_HAVE_GPU)
endif()

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
