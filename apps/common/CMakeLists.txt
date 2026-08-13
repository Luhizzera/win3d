if(WIN32)
    add_library(aether_app_common STATIC
        src/Gl33.cpp
        src/Ui.cpp
    )

    target_include_directories(aether_app_common PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )

    target_link_libraries(aether_app_common PUBLIC opengl32 gdi32 user32)
endif()
