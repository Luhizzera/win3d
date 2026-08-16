add_library(aether_solver STATIC
    src/DiffusionProblem.cpp
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
)

target_include_directories(aether_solver PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(aether_solver PUBLIC aether_mesh)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
