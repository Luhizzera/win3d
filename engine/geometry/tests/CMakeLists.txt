add_executable(aether_geometry_tests geometry_tests.cpp)
target_link_libraries(aether_geometry_tests PRIVATE aether_geometry aether_testing)

# The curved-solid round-trip test builds its own reference geometry with
# OpenCASCADE directly (BRepPrimAPI_MakeCylinder, STEPControl_Writer) to
# get a real curved STEP file without hand-authoring Part 21 text -- so
# this executable, not just aether_geometry, needs the same libraries and
# AETHER_HAVE_OPENCASCADE compile definition aether_geometry's own
# CMakeLists.txt sets PRIVATE for itself (which does not propagate here).
if(OpenCASCADE_FOUND)
    target_link_libraries(aether_geometry_tests PRIVATE
        TKDESTEP TKDEIGES TKXSBase TKDE
        TKMesh TKTopAlgo TKPrim TKBRep TKGeomBase TKG3d TKMath TKernel
    )
    target_compile_definitions(aether_geometry_tests PRIVATE AETHER_HAVE_OPENCASCADE)
endif()

add_test(NAME aether_geometry_tests COMMAND aether_geometry_tests)
