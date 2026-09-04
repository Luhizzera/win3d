add_library(aether_geometry STATIC
    src/TriangleMesh.cpp
    src/StlIO.cpp
    src/ObjIO.cpp
    src/StepIO.cpp
)

target_include_directories(aether_geometry PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# aether_mesh (not just aether_core) because StepIO reuses
# PolygonTriangulation2D to triangulate each planar STEP face -- no
# circular dependency, aether_mesh only depends on aether_core itself.
target_link_libraries(aether_geometry PUBLIC aether_core aether_mesh)

# Fase 5: curved-surface STEP/IGES, detected rather than required -- see
# the root CMakeLists.txt's own AETHER_ENABLE_OPENCASCADE comment. PRIVATE
# because AETHER_HAVE_OPENCASCADE only ever appears inside StepIO.cpp
# (never in a header StepIO.hpp's public callers see), the same rule
# AETHER_HAVE_GPU already follows for engine/solver -- no translation unit
# outside this library can ever disagree about whether the macro is
# defined. Library list is OpenCASCADE_DataExchange_LIBRARIES plus what it
# transitively needs for the tessellation/primitive-construction path
# (TKMesh, TKPrim) that DataExchange alone does not pull in.
if(OpenCASCADE_FOUND)
    target_link_libraries(aether_geometry PRIVATE
        TKDESTEP TKDEIGES TKXSBase TKDE
        TKMesh TKTopAlgo TKPrim TKBRep TKGeomBase TKG3d TKMath TKernel
    )
    target_compile_definitions(aether_geometry PRIVATE AETHER_HAVE_OPENCASCADE)
endif()

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
