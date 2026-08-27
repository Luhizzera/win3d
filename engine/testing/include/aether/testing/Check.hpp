#pragma once

#include <cstdio>
#include <cstdlib>

// A test-assertion macro that is always active, unlike <cassert>'s assert(),
// which becomes a no-op whenever NDEBUG is defined -- which CMake's default
// Release configuration does. Using plain assert() in test executables
// would make every "test" silently pass in Release builds regardless of
// what it actually checked.
#define AETHER_CHECK(condition)                                                                        \
    do {                                                                                                \
        if (!(condition)) {                                                                            \
            std::fprintf(stderr, "AETHER_CHECK failed: %s\n  at %s:%d\n", #condition, __FILE__, __LINE__); \
            std::abort();                                                                               \
        }                                                                                               \
    } while (false)
