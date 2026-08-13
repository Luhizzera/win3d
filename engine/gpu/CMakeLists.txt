add_library(aether_gpu STATIC
    src/PoissonOperatorCuda.cu
)

target_include_directories(aether_gpu PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Blackwell (this project's own dev machine, an RTX 5080) is compute
# capability 12.x, new enough that hardcoding an architecture list would
# need updating on the next GPU generation anyway -- "native" (needs CMake
# 3.24+; the root CMakeLists.txt's own cmake_minimum_required is lower, but
# every machine this has actually been built on so far has a newer CMake)
# asks nvcc to target whatever GPU is actually present at configure time.
set_target_properties(aether_gpu PROPERTIES CUDA_ARCHITECTURES native)

# nvcc fuses a*b+c into a single fused-multiply-add instruction by default
# (one rounding step instead of two) -- IEEE-754-legal, but NOT what plain
# scalar double-precision code on this project's MSVC target does without
# /fp:fast, so it silently breaks the bit-exact GPU-vs-CPU comparison this
# module's own tests rely on (measured directly: ~1e-13 to 1e-14 absolute
# differences on values up to a few hundred, before this flag -- exactly
# the signature of one extra/missing rounding, not an algorithmic bug).
# -fmad=false forces the same sequential-rounding arithmetic the CPU does.
target_compile_options(aether_gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-fmad=false>)

if(AETHER_BUILD_TESTS)
    add_subdirectory(tests)
endif()
