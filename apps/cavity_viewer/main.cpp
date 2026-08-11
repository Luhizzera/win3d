// Aether Cavity Viewer - Module 7 (post-processing): solves the lid-driven
// cavity (Module 4's LidDrivenCavitySolver2D) and renders both the raw
// velocity field as arrows *and* real integrated streamlines (RK4 via
// Streamline2D, not just arrows at cell centers) -- a direct visual check
// of the primary recirculating vortex the solver_tests.cpp topology test
// proves must exist (top row dragged by the lid, bottom row forced to
// reverse by mass conservation), now traced as actual particle paths.
// Modernized (third of the four legacy viewer apps) onto the shared
// OpenGL 3.3 core-profile / shader / VBO pipeline in apps/common/Gl33.

#include "aether/core/Matrix4x4.hpp"
#include "aether/postprocessing/Streamline2D.hpp"
#include "aether/solver/LidDrivenCavitySolver2D.hpp"
#include "aether_app/Gl33.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using aether::core::Matrix4x4;
using aether::solver::LidDrivenCavitySolver2D;
using aether::postprocessing::Streamline2D;
namespace gl33 = aether::app;

namespace {

// A single shader handles both the flat-colored arrow/outline geometry
// (GL_LINES, per-vertex color) and the streamline overlay (drawn in a
// second pass with a uniform color override folded in via vertex color).
const char* kVertexShaderSource = R"GLSL(
#version 330 core
layout(location = 0) in vec2 aPos;
layout(location = 1) in vec3 aColor;
uniform mat4 uProjection;
out vec3 vColor;
void main() {
    vColor = aColor;
    gl_Position = uProjection * vec4(aPos, 0.0, 1.0);
}
)GLSL";

const char* kFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
}
)GLSL";

int g_width = 800;
int g_height = 800;

// Blue (slow) -> white (mid) -> red (fast), for t in [0, 1].
void colorForValue(double t, float& r, float& g, float& b) {
    t = std::clamp(t, 0.0, 1.0);
    if (t < 0.5) {
        const float s = static_cast<float>(t / 0.5);
        r = s;
        g = s;
        b = 1.0f;
    } else {
        const float s = static_cast<float>((t - 0.5) / 0.5);
        r = 1.0f;
        g = 1.0f - s;
        b = 1.0f - s;
    }
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    case WM_SIZE:
        g_width = LOWORD(lParam);
        g_height = HIWORD(lParam);
        if (wglGetCurrentContext()) {
            glViewport(0, 0, g_width, g_height);
        }
        return 0;
    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
}

void pushVertex(std::vector<float>& out, double x, double y, float r, float g, float b) {
    out.push_back(static_cast<float>(x));
    out.push_back(static_cast<float>(y));
    out.push_back(r);
    out.push_back(g);
    out.push_back(b);
}

} // namespace

int main() {
    const std::size_t n = 24;
    const double length = 1.0;
    const double viscosity = 0.1;
    const double lidVelocity = 1.0; // Re = lidVelocity*length/viscosity = 10, safely laminar

    LidDrivenCavitySolver2D solver(n, n, length, length, viscosity, lidVelocity);
    const double dt = 0.3 * solver.stableTimeStep();
    const int steps = 800;
    for (int s = 0; s < steps; ++s) {
        solver.step(dt);
    }

    std::printf("Aether Cavity Viewer - cavidade com tampa deslizante (Re=%.0f)\n",
                lidVelocity * length / viscosity);
    std::printf("  grade %zux%zu, passos: %d, divergencia max: %.4f\n", n, n, steps, solver.maxDivergence());

    double maxSpeed = 0.0;
    std::vector<double> uField(n * n);
    std::vector<double> vField(n * n);
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double u = solver.u(i, j);
            const double v = solver.v(i, j);
            uField[i + j * n] = u;
            vField[i + j * n] = v;
            maxSpeed = std::max(maxSpeed, std::sqrt(u * u + v * v));
        }
    }

    // Real streamline integration (RK4, non-periodic since the cavity has
    // solid walls) from a grid of seed points spread across the interior.
    Streamline2D tracer(n, n, length, length, uField, vField, /*periodic=*/false);
    const double cellSizeForSeeds = length / static_cast<double>(n);
    std::vector<std::vector<aether::core::Vector3>> streamlines;
    for (int sj = 1; sj <= 4; ++sj) {
        for (int si = 1; si <= 4; ++si) {
            const double seedX = (static_cast<double>(si) / 5.0) * length;
            const double seedY = (static_cast<double>(sj) / 5.0) * length;
            streamlines.push_back(tracer.trace(seedX, seedY, cellSizeForSeeds * 0.15, 400));
        }
    }

    HWND hwnd = gl33::createSimpleWindow(L"AetherCavityViewerWindowClass", L"Aether Cavity Viewer", g_width,
                                          g_height, WndProc);
    if (!hwnd) {
        std::fprintf(stderr, "erro ao criar janela\n");
        return 1;
    }

    HDC hdc = nullptr;
    HGLRC glContext = gl33::createGl33Context(hwnd, hdc);
    if (!glContext) {
        return 1;
    }

    const GLuint vertexShader = gl33::compileShader(gl33::kGlVertexShader, kVertexShaderSource);
    const GLuint fragmentShader = gl33::compileShader(gl33::kGlFragmentShader, kFragmentShaderSource);
    const GLuint program = gl33::linkProgram(vertexShader, fragmentShader);
    gl33::glDeleteShader(vertexShader);
    gl33::glDeleteShader(fragmentShader);
    const GLint projectionLoc = gl33::glGetUniformLocation(program, "uProjection");

    // --- Static geometry, built once: cavity outline, arrow field. ---
    const double cellSize = length / static_cast<double>(n);
    const double arrowScale = cellSize * 0.9 / std::max(maxSpeed, 1e-6);

    std::vector<float> lineVertexData;
    // Outline.
    {
        const float oc[3] = {0.5f, 0.5f, 0.55f};
        const double corners[5][2] = {
            {0.0, 0.0}, {length, 0.0}, {length, length}, {0.0, length}, {0.0, 0.0}};
        for (int c = 0; c + 1 < 5; ++c) {
            pushVertex(lineVertexData, corners[c][0], corners[c][1], oc[0], oc[1], oc[2]);
            pushVertex(lineVertexData, corners[c + 1][0], corners[c + 1][1], oc[0], oc[1], oc[2]);
        }
    }
    // Arrows (shaft + two-line arrowhead), colored by speed.
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            const double cx = (static_cast<double>(i) + 0.5) * cellSize;
            const double cy = (static_cast<double>(j) + 0.5) * cellSize;
            const double u = solver.u(i, j);
            const double v = solver.v(i, j);
            const double speed = std::sqrt(u * u + v * v);

            float r = 0.0f, g = 0.0f, b = 0.0f;
            colorForValue(speed / std::max(maxSpeed, 1e-6), r, g, b);

            const double tipX = cx + u * arrowScale;
            const double tipY = cy + v * arrowScale;
            pushVertex(lineVertexData, cx, cy, r, g, b);
            pushVertex(lineVertexData, tipX, tipY, r, g, b);

            if (speed > 1e-6) {
                const double dirX = u / speed;
                const double dirY = v / speed;
                const double headLen = cellSize * 0.25;
                const double leftX = tipX - headLen * (dirX * 0.866 - dirY * 0.5);
                const double leftY = tipY - headLen * (dirY * 0.866 + dirX * 0.5);
                const double rightX = tipX - headLen * (dirX * 0.866 + dirY * 0.5);
                const double rightY = tipY - headLen * (dirY * 0.866 - dirX * 0.5);
                pushVertex(lineVertexData, tipX, tipY, r, g, b);
                pushVertex(lineVertexData, leftX, leftY, r, g, b);
                pushVertex(lineVertexData, tipX, tipY, r, g, b);
                pushVertex(lineVertexData, rightX, rightY, r, g, b);
            }
        }
    }
    const auto arrowFieldVertexCount = static_cast<GLsizei>(lineVertexData.size() / 5);

    // Streamlines appended after the arrow field, as GL_LINE_STRIP runs
    // (one draw call per streamline, since strip boundaries can't be
    // expressed within a single GL_LINES buffer).
    std::vector<GLint> streamlineFirst;
    std::vector<GLsizei> streamlineCount;
    for (const auto& path : streamlines) {
        streamlineFirst.push_back(static_cast<GLint>(lineVertexData.size() / 5));
        for (const auto& p : path) {
            pushVertex(lineVertexData, p.x, p.y, 1.0f, 0.9f, 0.2f);
        }
        streamlineCount.push_back(static_cast<GLsizei>(path.size()));
    }

    GLuint vao = 0;
    GLuint vbo = 0;
    gl33::glGenVertexArrays(1, &vao);
    gl33::glBindVertexArray(vao);
    gl33::glGenBuffers(1, &vbo);
    gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo);
    gl33::glBufferData(gl33::kGlArrayBuffer,
                        static_cast<gl33::GLsizeiptr>(lineVertexData.size() * sizeof(float)),
                        lineVertexData.data(), gl33::kGlStaticDraw);
    const GLsizei stride = 5 * sizeof(float);
    gl33::glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    gl33::glEnableVertexAttribArray(0);
    gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(2 * sizeof(float)));
    gl33::glEnableVertexAttribArray(1);
    gl33::glBindVertexArray(0);

    glClearColor(0.05f, 0.05f, 0.08f, 1.0f);
    ShowWindow(hwnd, SW_SHOW);
    glViewport(0, 0, g_width, g_height);

    const double margin = 0.05;
    const Matrix4x4 projection = Matrix4x4::ortho(-margin, length + margin, -margin, length + margin, -1.0, 1.0);

    MSG msg{};
    bool running = true;
    while (running) {
        while (PeekMessageW(&msg, nullptr, 0, 0, PM_REMOVE)) {
            if (msg.message == WM_QUIT) {
                running = false;
            }
            TranslateMessage(&msg);
            DispatchMessageW(&msg);
        }
        if (!running) {
            break;
        }

        glClear(GL_COLOR_BUFFER_BIT);

        gl33::glUseProgram(program);
        gl33::glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, projection.data());
        gl33::glBindVertexArray(vao);
        glDrawArrays(GL_LINES, 0, arrowFieldVertexCount);
        for (std::size_t s = 0; s < streamlineFirst.size(); ++s) {
            glDrawArrays(GL_LINE_STRIP, streamlineFirst[s], streamlineCount[s]);
        }
        gl33::glBindVertexArray(0);

        SwapBuffers(hdc);
    }

    gl33::glDeleteBuffers(1, &vbo);
    gl33::glDeleteVertexArrays(1, &vao);
    gl33::glDeleteProgram(program);
    wglMakeCurrent(nullptr, nullptr);
    wglDeleteContext(glContext);
    ReleaseDC(hwnd, hdc);
    DestroyWindow(hwnd);
    return 0;
}
