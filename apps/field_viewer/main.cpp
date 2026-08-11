// Aether Field Viewer - Module 7: solves the Module 5 steady heat
// conduction problem on a small 2D grid and renders the resulting
// temperature field as a blue-white-red heatmap. Modernized (second of the
// four legacy viewer apps, after apps/viewer) onto the shared OpenGL 3.3
// core-profile / shader / VBO pipeline in apps/common/Gl33 -- see that
// header, and apps/viewer's own comment, for the full rationale.

#include "aether/core/Matrix4x4.hpp"
#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/solver/SteadyDiffusionSolver.hpp"
#include "aether_app/Gl33.hpp"

#include <algorithm>
#include <cstdio>
#include <vector>

using aether::core::Matrix4x4;
using aether::core::Vector3;
using aether::mesh::StructuredGrid3D;
using aether::solver::SteadyDiffusionSolver;
namespace gl33 = aether::app;

namespace {

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
int g_height = 400;

// Blue (cold) -> white (mid) -> red (hot), for t in [0, 1].
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

} // namespace

int main() {
    const std::size_t nx = 40;
    const std::size_t ny = 20;
    const double lx = 2.0;
    const double ly = 1.0;
    const double tMin = 0.0;
    const double tMax = 100.0;

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(lx, ly, 0.1), nx, ny, 1);
    SteadyDiffusionSolver solver(grid);
    // Three sides at tMin, one (the top) at tMax: the classic 2D Laplace
    // "plate with one hot edge" problem.
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, tMin);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, tMin);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::YMin, tMin);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::YMax, tMax);
    const std::size_t iterations = solver.solveConjugateGradient();

    std::printf("Aether Field Viewer - conducao de calor estacionaria (Modulo 5) em grade %zux%zu\n", nx,
                ny);
    std::printf("  iteracoes ate convergir: %zu\n", iterations);

    HWND hwnd = gl33::createSimpleWindow(L"AetherFieldViewerWindowClass", L"Aether Field Viewer", g_width,
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

    // Each cell is 2 triangles (6 vertices), each vertex [x, y, r, g, b].
    // Built once: the field is solved once and never changes afterward, so
    // there is no need to rebuild this every frame.
    const Vector3 spacing = grid.spacing();
    std::vector<float> vertexData;
    vertexData.reserve(nx * ny * 6 * 5);
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const double t = solver.value(i, j, 0);
            const double normalized = (t - tMin) / (tMax - tMin);
            float r = 0.0f, g = 0.0f, b = 0.0f;
            colorForValue(normalized, r, g, b);

            const double x0 = i * spacing.x;
            const double x1 = (i + 1) * spacing.x;
            const double y0 = j * spacing.y;
            const double y1 = (j + 1) * spacing.y;

            const double quad[6][2] = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y0}, {x1, y1}, {x0, y1}};
            for (const auto& p : quad) {
                vertexData.push_back(static_cast<float>(p[0]));
                vertexData.push_back(static_cast<float>(p[1]));
                vertexData.push_back(r);
                vertexData.push_back(g);
                vertexData.push_back(b);
            }
        }
    }
    const auto vertexCount = static_cast<GLsizei>(nx * ny * 6);

    GLuint vao = 0;
    GLuint vbo = 0;
    gl33::glGenVertexArrays(1, &vao);
    gl33::glBindVertexArray(vao);
    gl33::glGenBuffers(1, &vbo);
    gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo);
    gl33::glBufferData(gl33::kGlArrayBuffer, static_cast<gl33::GLsizeiptr>(vertexData.size() * sizeof(float)),
                        vertexData.data(), gl33::kGlStaticDraw);
    const GLsizei stride = 5 * sizeof(float);
    gl33::glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    gl33::glEnableVertexAttribArray(0);
    gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(2 * sizeof(float)));
    gl33::glEnableVertexAttribArray(1);
    gl33::glBindVertexArray(0);

    glClearColor(0.05f, 0.05f, 0.08f, 1.0f);

    ShowWindow(hwnd, SW_SHOW);
    glViewport(0, 0, g_width, g_height);

    const Matrix4x4 projection = Matrix4x4::ortho(0.0, lx, 0.0, ly, -1.0, 1.0);

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
        glDrawArrays(GL_TRIANGLES, 0, vertexCount);
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
