// Aether Unified Viewer - Module 8: a single executable combining the four
// previously-separate GL 3.3 core-profile viewer apps (viewer, field_viewer,
// cavity_viewer, turbulence_viewer), selected by a command-line mode
// argument. This is the literal "unificacao dos viewers" clause of the
// roadmap's modern-visualizer item, left open when each app was migrated
// individually onto the shared apps/common/Gl33 toolkit.
//
// Deliberately *not* a deeper architectural merge: each mode keeps its own
// self-contained setup/shader/VBO/message-loop code (wrapped in its own
// namespace below, copied near-verbatim from the four original apps) rather
// than forcing a single shared render-loop abstraction across genuinely
// different rendering needs (a 3D perspective scene with an interactive
// orbit camera vs. three static 2D orthographic overlays). That would add
// real complexity for little benefit, since only one mode ever runs per
// process invocation anyway -- this is "one executable, mode-selected", not
// "one code path for everything". The four original standalone apps were
// later removed in the cleanup before modules 9-14: keeping them meant
// carrying a second copy of the same render code, bound to drift.
//
// Usage:
//   aether_unified_viewer.exe mesh <arquivo.stl>
//   aether_unified_viewer.exe heatmap
//   aether_unified_viewer.exe cavity
//   aether_unified_viewer.exe cavity3d
//   aether_unified_viewer.exe isosurface
//   aether_unified_viewer.exe turbulence
//
// cavity3d (added later, closing a real Module-8 gap): every 3D solver
// built after this file's original four modes shipped (StaggeredLidDrivenCavitySolver3D
// and its three turbulence-closure extensions) had no viewer at all --
// console diagnostics only. cavity3d renders StaggeredLidDrivenCavitySolver3D
// itself (the shared base every 3D closure extends) with a 3D velocity
// vector field, reusing the mesh mode's orbit-camera/perspective skeleton.
// Deliberately does not add a mode per 3D closure (mixing-length/k-epsilon/
// k-omega SST 3D) -- they'd all render identically (same momentum solver,
// same vector field shape) and the interesting *closure-specific* fields
// (nu_t, k, epsilon/omega) are scalar volumes better suited to an
// iso-surface/slice renderer than arrows, which needs marching cubes 3D
// (not yet built -- engine/postprocessing only has the 2D marching squares
// this project's own docs already flag as "the first step toward" it) to
// be worth building well rather than half-built here.
//
// isosurface (added once marching cubes 3D existed): renders the eddy-
// viscosity (nu_t) field of KOmegaSSTLidDrivenCavitySolver3D -- this
// project's most complete 3D turbulence closure -- as a real shaded
// iso-surface via aether::postprocessing::marchingCubes3D, reusing the
// mesh mode's flat-shaded lit-triangle pipeline (position+normal, two-
// sided lambertian) rather than the arrow-based cavity3d mode, since a
// scalar volume is what marching cubes was built for.

#include "aether/core/Matrix4x4.hpp"
#include "aether/geometry/StlIO.hpp"
#include "aether/geometry/TriangleMesh.hpp"
#include "aether/mesh/StructuredGrid3D.hpp"
#include "aether/postprocessing/MarchingCubes3D.hpp"
#include "aether/postprocessing/Streamline2D.hpp"
#include "aether/solver/KEpsilonChannelFlowSolver1D.hpp"
#include "aether/solver/KOmegaSSTChannelFlowSolver1D.hpp"
#include "aether/solver/KOmegaSSTLidDrivenCavitySolver3D.hpp"
#include "aether/solver/KEpsilonLidDrivenCavitySolver2D.hpp"
#include "aether/solver/KOmegaSSTLidDrivenCavitySolver2D.hpp"
#include "aether/solver/LidDrivenCavitySolver2D.hpp"
#include "aether/solver/MixingLengthLidDrivenCavitySolver2D.hpp"
#include "aether/solver/MixingLengthChannelFlowSolver1D.hpp"
#include "aether/solver/StaggeredLidDrivenCavitySolver3D.hpp"
#include "aether/solver/SteadyDiffusionSolver.hpp"
#include "aether_app/Gl33.hpp"
#include "aether_app/Ui.hpp"

#include <windowsx.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <memory>
#include <string>
#include <vector>

namespace gl33 = aether::app;

// ===========================================================================
// Shared helpers used by more than one mode.
// ===========================================================================
namespace {

// Blue (low) -> white (mid) -> red (high), for t in [0, 1]. Used by both the
// heatmap and cavity modes for their respective color scales.
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

void pushVertex(std::vector<float>& out, double x, double y, float r, float g, float b) {
    out.push_back(static_cast<float>(x));
    out.push_back(static_cast<float>(y));
    out.push_back(r);
    out.push_back(g);
    out.push_back(b);
}

} // namespace

// ===========================================================================
// Mode: mesh (formerly apps/viewer) -- STL mesh with an orbit camera.
// ===========================================================================
namespace mesh_mode {

using aether::core::Matrix4x4;
using aether::core::Vector3;
using aether::geometry::loadStl;
using aether::geometry::TriangleMesh;

const char* kVertexShaderSource = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
uniform mat4 uMvp;
out vec3 vNormal;
void main() {
    vNormal = aNormal;
    gl_Position = uMvp * vec4(aPos, 1.0);
}
)GLSL";

const char* kFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vNormal;
uniform vec3 uLightDir;
uniform vec3 uBaseColor;
out vec4 FragColor;
void main() {
    vec3 n = normalize(vNormal);
    float diff = abs(dot(n, normalize(uLightDir)));
    vec3 color = uBaseColor * (0.35 + 0.65 * diff);
    FragColor = vec4(color, 1.0);
}
)GLSL";

struct OrbitCamera {
    Vector3 center;
    double distance = 5.0;
    double yaw = 0.7;
    double pitch = 0.5;

    Vector3 eye() const {
        const double cp = std::cos(pitch);
        return center + Vector3(distance * cp * std::cos(yaw), distance * std::sin(pitch),
                                 distance * cp * std::sin(yaw));
    }
};

OrbitCamera g_camera;
bool g_dragging = false;
POINT g_lastMouse{};
int g_width = 1024;
int g_height = 768;

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
    case WM_LBUTTONDOWN:
        g_dragging = true;
        g_lastMouse.x = GET_X_LPARAM(lParam);
        g_lastMouse.y = GET_Y_LPARAM(lParam);
        SetCapture(hwnd);
        return 0;
    case WM_LBUTTONUP:
        g_dragging = false;
        ReleaseCapture();
        return 0;
    case WM_MOUSEMOVE:
        if (g_dragging) {
            const int x = GET_X_LPARAM(lParam);
            const int y = GET_Y_LPARAM(lParam);
            g_camera.yaw += (x - g_lastMouse.x) * 0.01;
            g_camera.pitch += (y - g_lastMouse.y) * 0.01;
            g_camera.pitch = std::clamp(g_camera.pitch, -1.5, 1.5);
            g_lastMouse.x = x;
            g_lastMouse.y = y;
        }
        return 0;
    case WM_MOUSEWHEEL: {
        const short delta = GET_WHEEL_DELTA_WPARAM(wParam);
        const double factor = delta > 0 ? 0.9 : 1.1;
        g_camera.distance = std::max(1e-3, g_camera.distance * factor);
        return 0;
    }
    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
}

int run(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "uso: aether_unified_viewer mesh <arquivo.stl>\n");
        return 1;
    }

    TriangleMesh mesh;
    try {
        mesh = loadStl(argv[2]);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "erro ao carregar STL: %s\n", e.what());
        return 1;
    }

    Vector3 minBound = mesh.vertexCount() > 0 ? mesh.vertex(0) : Vector3{};
    Vector3 maxBound = minBound;
    for (std::size_t i = 1; i < mesh.vertexCount(); ++i) {
        const Vector3& p = mesh.vertex(i);
        minBound.x = std::min(minBound.x, p.x);
        maxBound.x = std::max(maxBound.x, p.x);
        minBound.y = std::min(minBound.y, p.y);
        maxBound.y = std::max(maxBound.y, p.y);
        minBound.z = std::min(minBound.z, p.z);
        maxBound.z = std::max(maxBound.z, p.z);
    }
    g_camera.center = (minBound + maxBound) * 0.5;
    const Vector3 extent = maxBound - minBound;
    const double diagonal = std::sqrt(extent.x * extent.x + extent.y * extent.y + extent.z * extent.z);
    g_camera.distance = diagonal > 0 ? diagonal * 1.5 : 5.0;

    std::printf("Aether Unified Viewer [mesh] - %s\n", argv[2]);
    std::printf("  vertices: %zu, triangulos: %zu\n", mesh.vertexCount(), mesh.triangleCount());
    std::printf("  area: %.6f, volume: %.6f, watertight: %s\n", mesh.surfaceArea(), mesh.volume(),
                mesh.isWatertight() ? "sim" : "nao");

    HWND hwnd = gl33::createSimpleWindow(L"AetherUnifiedViewerMeshWindowClass", L"Aether Viewer - Mesh",
                                          g_width, g_height, WndProc);
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

    const GLint mvpLoc = gl33::glGetUniformLocation(program, "uMvp");
    const GLint lightDirLoc = gl33::glGetUniformLocation(program, "uLightDir");
    const GLint baseColorLoc = gl33::glGetUniformLocation(program, "uBaseColor");

    std::vector<float> vertexData;
    vertexData.reserve(mesh.triangleCount() * 3 * 6);
    for (std::size_t t = 0; t < mesh.triangleCount(); ++t) {
        const Vector3 n = mesh.faceNormal(t);
        const auto& tri = mesh.triangle(t);
        for (std::size_t v = 0; v < 3; ++v) {
            const Vector3& p = mesh.vertex(tri.vertices[v]);
            vertexData.push_back(static_cast<float>(p.x));
            vertexData.push_back(static_cast<float>(p.y));
            vertexData.push_back(static_cast<float>(p.z));
            vertexData.push_back(static_cast<float>(n.x));
            vertexData.push_back(static_cast<float>(n.y));
            vertexData.push_back(static_cast<float>(n.z));
        }
    }
    const auto vertexCount = static_cast<GLsizei>(mesh.triangleCount() * 3);

    GLuint vao = 0;
    GLuint vbo = 0;
    gl33::glGenVertexArrays(1, &vao);
    gl33::glBindVertexArray(vao);
    gl33::glGenBuffers(1, &vbo);
    gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo);
    gl33::glBufferData(gl33::kGlArrayBuffer, static_cast<gl33::GLsizeiptr>(vertexData.size() * sizeof(float)),
                        vertexData.data(), gl33::kGlStaticDraw);
    const GLsizei stride = 6 * sizeof(float);
    gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    gl33::glEnableVertexAttribArray(0);
    gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));
    gl33::glEnableVertexAttribArray(1);
    gl33::glBindVertexArray(0);

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.1f, 0.1f, 0.15f, 1.0f);

    ShowWindow(hwnd, SW_SHOW);
    glViewport(0, 0, g_width, g_height);

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

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const double aspect = g_height > 0 ? static_cast<double>(g_width) / g_height : 1.0;
        const double nearPlane = std::max(0.01 * g_camera.distance, 1e-3);
        const double farPlane = g_camera.distance * 100.0 + 10.0;
        const Matrix4x4 projection = Matrix4x4::perspective(45.0 * 3.14159265358979323846 / 180.0, aspect,
                                                              nearPlane, farPlane);
        const Vector3 eye = g_camera.eye();
        const Matrix4x4 view = Matrix4x4::lookAt(eye, g_camera.center, Vector3(0.0, 1.0, 0.0));
        const Matrix4x4 mvp = projection * view;

        gl33::glUseProgram(program);
        gl33::glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, mvp.data());
        gl33::glUniform3f(lightDirLoc, 1.0f, 1.0f, 1.0f);
        gl33::glUniform3f(baseColorLoc, 0.7f, 0.75f, 0.85f);

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

} // namespace mesh_mode

// ===========================================================================
// Mode: heatmap (formerly apps/field_viewer) -- steady heat conduction.
// ===========================================================================
namespace heatmap_mode {

using aether::core::Matrix4x4;
using aether::core::Vector3;
using aether::mesh::StructuredGrid3D;
using aether::solver::SteadyDiffusionSolver;

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

int run() {
    const std::size_t nx = 40;
    const std::size_t ny = 20;
    const double lx = 2.0;
    const double ly = 1.0;
    const double tMin = 0.0;
    const double tMax = 100.0;

    StructuredGrid3D grid(Vector3(0.0, 0.0, 0.0), Vector3(lx, ly, 0.1), nx, ny, 1);
    SteadyDiffusionSolver solver(grid);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMin, tMin);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::XMax, tMin);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::YMin, tMin);
    solver.setBoundaryValue(SteadyDiffusionSolver::Face::YMax, tMax);
    const std::size_t iterations = solver.solveConjugateGradient();

    std::printf("Aether Unified Viewer [heatmap] - conducao de calor estacionaria em grade %zux%zu\n", nx,
                ny);
    std::printf("  iteracoes ate convergir: %zu\n", iterations);

    HWND hwnd = gl33::createSimpleWindow(L"AetherUnifiedViewerHeatmapWindowClass",
                                          L"Aether Viewer - Heatmap", g_width, g_height, WndProc);
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

} // namespace heatmap_mode

// ===========================================================================
// Mode: cavity (formerly apps/cavity_viewer) -- velocity arrows + streamlines.
// ===========================================================================
namespace cavity_mode {

using aether::core::Matrix4x4;
using aether::solver::LidDrivenCavitySolver2D;
using aether::postprocessing::Streamline2D;

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

int run() {
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

    std::printf("Aether Unified Viewer [cavity] - cavidade com tampa deslizante (Re=%.0f)\n",
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

    HWND hwnd = gl33::createSimpleWindow(L"AetherUnifiedViewerCavityWindowClass", L"Aether Viewer - Cavity",
                                          g_width, g_height, WndProc);
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

    const double cellSize = length / static_cast<double>(n);
    const double arrowScale = cellSize * 0.9 / std::max(maxSpeed, 1e-6);

    std::vector<float> lineVertexData;
    {
        const float oc[3] = {0.5f, 0.5f, 0.55f};
        const double corners[5][2] = {
            {0.0, 0.0}, {length, 0.0}, {length, length}, {0.0, length}, {0.0, 0.0}};
        for (int c = 0; c + 1 < 5; ++c) {
            pushVertex(lineVertexData, corners[c][0], corners[c][1], oc[0], oc[1], oc[2]);
            pushVertex(lineVertexData, corners[c + 1][0], corners[c + 1][1], oc[0], oc[1], oc[2]);
        }
    }
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

} // namespace cavity_mode

// ===========================================================================
// Mode: cavity3d -- the first viewer for any of this project's 3D solvers.
// StaggeredLidDrivenCavitySolver3D's velocity field as a 3D arrow field,
// colored by speed, on the mesh mode's perspective/orbit-camera skeleton.
// ===========================================================================
namespace cavity3d_mode {

using aether::core::Matrix4x4;
using aether::core::Vector3;
using aether::solver::StaggeredLidDrivenCavitySolver3D;

const char* kVertexShaderSource = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aColor;
uniform mat4 uMvp;
out vec3 vColor;
void main() {
    vColor = aColor;
    gl_Position = uMvp * vec4(aPos, 1.0);
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

struct OrbitCamera {
    Vector3 center;
    double distance = 3.0;
    double yaw = 0.7;
    double pitch = 0.5;

    Vector3 eye() const {
        const double cp = std::cos(pitch);
        return center + Vector3(distance * cp * std::cos(yaw), distance * std::sin(pitch),
                                 distance * cp * std::sin(yaw));
    }
};

OrbitCamera g_camera;
bool g_dragging = false;
POINT g_lastMouse{};
int g_width = 1024;
int g_height = 768;

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
    case WM_LBUTTONDOWN:
        g_dragging = true;
        g_lastMouse.x = GET_X_LPARAM(lParam);
        g_lastMouse.y = GET_Y_LPARAM(lParam);
        SetCapture(hwnd);
        return 0;
    case WM_LBUTTONUP:
        g_dragging = false;
        ReleaseCapture();
        return 0;
    case WM_MOUSEMOVE:
        if (g_dragging) {
            const int x = GET_X_LPARAM(lParam);
            const int y = GET_Y_LPARAM(lParam);
            g_camera.yaw += (x - g_lastMouse.x) * 0.01;
            g_camera.pitch += (y - g_lastMouse.y) * 0.01;
            g_camera.pitch = std::clamp(g_camera.pitch, -1.5, 1.5);
            g_lastMouse.x = x;
            g_lastMouse.y = y;
        }
        return 0;
    case WM_MOUSEWHEEL: {
        const short delta = GET_WHEEL_DELTA_WPARAM(wParam);
        const double factor = delta > 0 ? 0.9 : 1.1;
        g_camera.distance = std::max(1e-3, g_camera.distance * factor);
        return 0;
    }
    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
}

int run() {
    const std::size_t n = 10;
    const double length = 1.0;
    const double viscosity = 0.1;
    const double lidVelocity = 1.0; // Re = lidVelocity*length/viscosity = 10, matches this
                                     // class's own validated test case

    StaggeredLidDrivenCavitySolver3D solver(n, n, n, length, length, length, viscosity, lidVelocity);
    const double dt = 0.3 * solver.stableTimeStep();
    const int steps = 400;
    for (int s = 0; s < steps; ++s) {
        solver.step(dt);
    }

    std::printf("Aether Unified Viewer [cavity3d] - cavidade 3D com tampa deslizante (Re=%.0f)\n",
                lidVelocity * length / viscosity);
    std::printf("  grade %zux%zux%zu, passos: %d, divergencia max: %.6f\n", n, n, n, steps,
                solver.maxDivergence());

    g_camera.center = Vector3(length * 0.5, length * 0.5, length * 0.5);
    g_camera.distance = length * 2.5;

    HWND hwnd = gl33::createSimpleWindow(L"AetherUnifiedViewerCavity3DWindowClass",
                                          L"Aether Viewer - Cavity 3D", g_width, g_height, WndProc);
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
    const GLint mvpLoc = gl33::glGetUniformLocation(program, "uMvp");

    // Cell-centered velocity, interpolated from the staggered faces
    // (the same 0.5*(uAt(i)+uAt(i+1)) style average every 3D solver in
    // this project already uses to sample velocity at a cell center).
    double maxSpeed = 0.0;
    std::vector<double> speedAt(n * n * n);
    auto flatIndex = [n](std::size_t i, std::size_t j, std::size_t k) { return i + j * n + k * n * n; };
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double uc = 0.5 * (solver.u(i, j, k) + solver.u(i + 1, j, k));
                const double vc = 0.5 * (solver.v(i, j, k) + solver.v(i, j + 1, k));
                const double wc = 0.5 * (solver.w(i, j, k) + solver.w(i, j, k + 1));
                const double speed = std::sqrt(uc * uc + vc * vc + wc * wc);
                speedAt[flatIndex(i, j, k)] = speed;
                maxSpeed = std::max(maxSpeed, speed);
            }
        }
    }

    const double cellSize = length / static_cast<double>(n);
    const double arrowScale = cellSize * 0.9 / std::max(maxSpeed, 1e-6);

    std::vector<float> vertexData;
    auto pushVertex3 = [&vertexData](double x, double y, double z, float r, float g, float b) {
        vertexData.push_back(static_cast<float>(x));
        vertexData.push_back(static_cast<float>(y));
        vertexData.push_back(static_cast<float>(z));
        vertexData.push_back(r);
        vertexData.push_back(g);
        vertexData.push_back(b);
    };

    // Wireframe box, for spatial reference.
    {
        const float c[3] = {0.5f, 0.5f, 0.55f};
        const double corners[8][3] = {
            {0, 0, 0}, {length, 0, 0}, {length, length, 0}, {0, length, 0},
            {0, 0, length}, {length, 0, length}, {length, length, length}, {0, length, length},
        };
        const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                   {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        for (const auto& e : edges) {
            pushVertex3(corners[e[0]][0], corners[e[0]][1], corners[e[0]][2], c[0], c[1], c[2]);
            pushVertex3(corners[e[1]][0], corners[e[1]][1], corners[e[1]][2], c[0], c[1], c[2]);
        }
    }
    const auto boxVertexCount = static_cast<GLsizei>(vertexData.size() / 6);

    // Velocity arrows (shafts only -- a 3D arrowhead needs an orientation
    // frame perpendicular to the shaft, real extra complexity for little
    // added readability once color-by-speed and line direction already
    // convey the field; the 2D cavity mode's arrowheads work cheaply
    // because "perpendicular to a 2D line" has only one choice, which
    // doesn't generalize to 3D).
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double cx = (static_cast<double>(i) + 0.5) * cellSize;
                const double cy = (static_cast<double>(j) + 0.5) * cellSize;
                const double cz = (static_cast<double>(k) + 0.5) * cellSize;
                const double uc = 0.5 * (solver.u(i, j, k) + solver.u(i + 1, j, k));
                const double vc = 0.5 * (solver.v(i, j, k) + solver.v(i, j + 1, k));
                const double wc = 0.5 * (solver.w(i, j, k) + solver.w(i, j, k + 1));
                const double speed = speedAt[flatIndex(i, j, k)];

                float r = 0.0f, g = 0.0f, b = 0.0f;
                colorForValue(speed / std::max(maxSpeed, 1e-6), r, g, b);

                pushVertex3(cx, cy, cz, r, g, b);
                pushVertex3(cx + uc * arrowScale, cy + vc * arrowScale, cz + wc * arrowScale, r, g, b);
            }
        }
    }
    const auto arrowVertexCount = static_cast<GLsizei>(vertexData.size() / 6 - boxVertexCount);

    GLuint vao = 0;
    GLuint vbo = 0;
    gl33::glGenVertexArrays(1, &vao);
    gl33::glBindVertexArray(vao);
    gl33::glGenBuffers(1, &vbo);
    gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo);
    gl33::glBufferData(gl33::kGlArrayBuffer, static_cast<gl33::GLsizeiptr>(vertexData.size() * sizeof(float)),
                        vertexData.data(), gl33::kGlStaticDraw);
    const GLsizei stride = 6 * sizeof(float);
    gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    gl33::glEnableVertexAttribArray(0);
    gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));
    gl33::glEnableVertexAttribArray(1);
    gl33::glBindVertexArray(0);

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.08f, 0.08f, 0.1f, 1.0f);
    glLineWidth(1.5f);

    ShowWindow(hwnd, SW_SHOW);
    glViewport(0, 0, g_width, g_height);

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

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const double aspect = g_height > 0 ? static_cast<double>(g_width) / g_height : 1.0;
        const double nearPlane = std::max(0.01 * g_camera.distance, 1e-3);
        const double farPlane = g_camera.distance * 100.0 + 10.0;
        const Matrix4x4 projection = Matrix4x4::perspective(45.0 * 3.14159265358979323846 / 180.0, aspect,
                                                              nearPlane, farPlane);
        const Vector3 eye = g_camera.eye();
        const Matrix4x4 view = Matrix4x4::lookAt(eye, g_camera.center, Vector3(0.0, 1.0, 0.0));
        const Matrix4x4 mvp = projection * view;

        gl33::glUseProgram(program);
        gl33::glUniformMatrix4fv(mvpLoc, 1, GL_FALSE, mvp.data());

        gl33::glBindVertexArray(vao);
        glDrawArrays(GL_LINES, 0, boxVertexCount);
        glDrawArrays(GL_LINES, boxVertexCount, arrowVertexCount);
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

} // namespace cavity3d_mode

// ===========================================================================
// Mode: isosurface -- eddy-viscosity iso-surface from the 3D k-omega SST
// cavity, extracted via marchingCubes3D and rendered as a real shaded
// triangle mesh (mesh mode's lit pipeline), not arrows.
// ===========================================================================
namespace isosurface_mode {

using aether::core::Matrix4x4;
using aether::core::Vector3;
using aether::postprocessing::marchingCubes3D;
using aether::postprocessing::Triangle3D;
using aether::solver::KOmegaSSTLidDrivenCavitySolver3D;

const char* kSurfaceVertexShaderSource = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNormal;
uniform mat4 uMvp;
out vec3 vNormal;
void main() {
    vNormal = aNormal;
    gl_Position = uMvp * vec4(aPos, 1.0);
}
)GLSL";

const char* kSurfaceFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vNormal;
uniform vec3 uLightDir;
uniform vec3 uBaseColor;
out vec4 FragColor;
void main() {
    vec3 n = normalize(vNormal);
    float diff = abs(dot(n, normalize(uLightDir)));
    vec3 color = uBaseColor * (0.35 + 0.65 * diff);
    FragColor = vec4(color, 1.0);
}
)GLSL";

const char* kLineVertexShaderSource = R"GLSL(
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aColor;
uniform mat4 uMvp;
out vec3 vColor;
void main() {
    vColor = aColor;
    gl_Position = uMvp * vec4(aPos, 1.0);
}
)GLSL";

const char* kLineFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
}
)GLSL";

struct OrbitCamera {
    Vector3 center;
    double distance = 3.0;
    double yaw = 0.7;
    double pitch = 0.5;

    Vector3 eye() const {
        const double cp = std::cos(pitch);
        return center + Vector3(distance * cp * std::cos(yaw), distance * std::sin(pitch),
                                 distance * cp * std::sin(yaw));
    }
};

OrbitCamera g_camera;
bool g_dragging = false;
POINT g_lastMouse{};
int g_width = 1024;
int g_height = 768;

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
    case WM_LBUTTONDOWN:
        g_dragging = true;
        g_lastMouse.x = GET_X_LPARAM(lParam);
        g_lastMouse.y = GET_Y_LPARAM(lParam);
        SetCapture(hwnd);
        return 0;
    case WM_LBUTTONUP:
        g_dragging = false;
        ReleaseCapture();
        return 0;
    case WM_MOUSEMOVE:
        if (g_dragging) {
            const int x = GET_X_LPARAM(lParam);
            const int y = GET_Y_LPARAM(lParam);
            g_camera.yaw += (x - g_lastMouse.x) * 0.01;
            g_camera.pitch += (y - g_lastMouse.y) * 0.01;
            g_camera.pitch = std::clamp(g_camera.pitch, -1.5, 1.5);
            g_lastMouse.x = x;
            g_lastMouse.y = y;
        }
        return 0;
    case WM_MOUSEWHEEL: {
        const short delta = GET_WHEEL_DELTA_WPARAM(wParam);
        const double factor = delta > 0 ? 0.9 : 1.1;
        g_camera.distance = std::max(1e-3, g_camera.distance * factor);
        return 0;
    }
    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
}

int run() {
    // n=12/300 steps deliberately matches (not exceeds) the resolution
    // KOmegaSSTLidDrivenCavitySolver3D's own C++ test already runs at --
    // this class's warm-start (400 internal mixing-length steps) plus its
    // own transport equations make it by far the most expensive solver any
    // viewer mode in this project drives; a Debug build measured ~8s total
    // at this size versus ~32s at n=16/500 steps, and an interactive
    // viewer should not make the user wait half a minute before the first
    // frame.
    const std::size_t n = 12;
    const double length = 1.0;
    const double viscosity = 0.01;
    const double lidVelocity = 1.0; // Re = 100, matching this class's own validated test case

    KOmegaSSTLidDrivenCavitySolver3D solver(n, n, n, length, length, length, viscosity, lidVelocity);
    double dt = 0.3 * solver.stableTimeStep();
    const int steps = 300;
    for (int s = 0; s < steps; ++s) {
        solver.step(dt);
        if (s % 100 == 0) {
            dt = 0.3 * solver.stableTimeStep();
        }
    }

    std::vector<double> nut(n * n * n);
    double maxNut = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const double v = solver.eddyViscosity(i, j, k);
                nut[i + j * n + k * n * n] = v;
                maxNut = std::max(maxNut, v);
            }
        }
    }

    // Relative threshold rather than a fixed absolute isoValue: nu_t's
    // magnitude is problem-dependent (viscosity, Re, resolution all shift
    // it), so a fraction of the observed maximum guarantees a non-empty,
    // reasonably-sized surface regardless of the exact numbers this run
    // produces.
    const double isoValue = 0.3 * std::max(maxNut, 1e-12);
    const auto triangles = marchingCubes3D(n, n, n, length, length, length, nut, isoValue);

    std::printf("Aether Unified Viewer [isosurface] - iso-superficie de nu_t (k-omega SST 3D, Re=%.0f)\n",
                lidVelocity * length / viscosity);
    std::printf("  grade %zux%zux%zu, passos: %d, max(nu_t)=%.6e, iso=%.6e, triangulos: %zu\n", n, n, n,
                steps, maxNut, isoValue, triangles.size());

    g_camera.center = Vector3(length * 0.5, length * 0.5, length * 0.5);
    g_camera.distance = length * 2.5;

    HWND hwnd = gl33::createSimpleWindow(L"AetherUnifiedViewerIsosurfaceWindowClass",
                                          L"Aether Viewer - Isosurface", g_width, g_height, WndProc);
    if (!hwnd) {
        std::fprintf(stderr, "erro ao criar janela\n");
        return 1;
    }

    HDC hdc = nullptr;
    HGLRC glContext = gl33::createGl33Context(hwnd, hdc);
    if (!glContext) {
        return 1;
    }

    const GLuint surfaceVs = gl33::compileShader(gl33::kGlVertexShader, kSurfaceVertexShaderSource);
    const GLuint surfaceFs = gl33::compileShader(gl33::kGlFragmentShader, kSurfaceFragmentShaderSource);
    const GLuint surfaceProgram = gl33::linkProgram(surfaceVs, surfaceFs);
    gl33::glDeleteShader(surfaceVs);
    gl33::glDeleteShader(surfaceFs);
    const GLint surfaceMvpLoc = gl33::glGetUniformLocation(surfaceProgram, "uMvp");
    const GLint lightDirLoc = gl33::glGetUniformLocation(surfaceProgram, "uLightDir");
    const GLint baseColorLoc = gl33::glGetUniformLocation(surfaceProgram, "uBaseColor");

    const GLuint lineVs = gl33::compileShader(gl33::kGlVertexShader, kLineVertexShaderSource);
    const GLuint lineFs = gl33::compileShader(gl33::kGlFragmentShader, kLineFragmentShaderSource);
    const GLuint lineProgram = gl33::linkProgram(lineVs, lineFs);
    gl33::glDeleteShader(lineVs);
    gl33::glDeleteShader(lineFs);
    const GLint lineMvpLoc = gl33::glGetUniformLocation(lineProgram, "uMvp");

    // Surface mesh: flat per-face shading, same convention as mesh_mode's
    // STL rendering (each triangle's 3 vertices repeat its own face
    // normal).
    std::vector<float> surfaceVertexData;
    surfaceVertexData.reserve(triangles.size() * 3 * 6);
    for (const Triangle3D& tri : triangles) {
        const Vector3 normal = (tri.b - tri.a).cross(tri.c - tri.a).normalized();
        for (const Vector3& p : {tri.a, tri.b, tri.c}) {
            surfaceVertexData.push_back(static_cast<float>(p.x));
            surfaceVertexData.push_back(static_cast<float>(p.y));
            surfaceVertexData.push_back(static_cast<float>(p.z));
            surfaceVertexData.push_back(static_cast<float>(normal.x));
            surfaceVertexData.push_back(static_cast<float>(normal.y));
            surfaceVertexData.push_back(static_cast<float>(normal.z));
        }
    }
    const auto surfaceVertexCount = static_cast<GLsizei>(triangles.size() * 3);

    GLuint surfaceVao = 0;
    GLuint surfaceVbo = 0;
    gl33::glGenVertexArrays(1, &surfaceVao);
    gl33::glBindVertexArray(surfaceVao);
    gl33::glGenBuffers(1, &surfaceVbo);
    gl33::glBindBuffer(gl33::kGlArrayBuffer, surfaceVbo);
    gl33::glBufferData(gl33::kGlArrayBuffer,
                        static_cast<gl33::GLsizeiptr>(surfaceVertexData.size() * sizeof(float)),
                        surfaceVertexData.empty() ? nullptr : surfaceVertexData.data(), gl33::kGlStaticDraw);
    const GLsizei surfaceStride = 6 * sizeof(float);
    gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, surfaceStride, reinterpret_cast<void*>(0));
    gl33::glEnableVertexAttribArray(0);
    gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, surfaceStride,
                                 reinterpret_cast<void*>(3 * sizeof(float)));
    gl33::glEnableVertexAttribArray(1);
    gl33::glBindVertexArray(0);

    // Wireframe box, for spatial reference -- same pattern as cavity3d_mode.
    std::vector<float> lineVertexData;
    {
        const float c[3] = {0.5f, 0.5f, 0.55f};
        const double corners[8][3] = {
            {0, 0, 0}, {length, 0, 0}, {length, length, 0}, {0, length, 0},
            {0, 0, length}, {length, 0, length}, {length, length, length}, {0, length, length},
        };
        const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                   {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        for (const auto& e : edges) {
            for (int v = 0; v < 2; ++v) {
                lineVertexData.push_back(static_cast<float>(corners[e[v]][0]));
                lineVertexData.push_back(static_cast<float>(corners[e[v]][1]));
                lineVertexData.push_back(static_cast<float>(corners[e[v]][2]));
                lineVertexData.push_back(c[0]);
                lineVertexData.push_back(c[1]);
                lineVertexData.push_back(c[2]);
            }
        }
    }
    const auto boxVertexCount = static_cast<GLsizei>(lineVertexData.size() / 6);

    GLuint lineVao = 0;
    GLuint lineVbo = 0;
    gl33::glGenVertexArrays(1, &lineVao);
    gl33::glBindVertexArray(lineVao);
    gl33::glGenBuffers(1, &lineVbo);
    gl33::glBindBuffer(gl33::kGlArrayBuffer, lineVbo);
    gl33::glBufferData(gl33::kGlArrayBuffer,
                        static_cast<gl33::GLsizeiptr>(lineVertexData.size() * sizeof(float)),
                        lineVertexData.data(), gl33::kGlStaticDraw);
    const GLsizei lineStride = 6 * sizeof(float);
    gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, lineStride, reinterpret_cast<void*>(0));
    gl33::glEnableVertexAttribArray(0);
    gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, lineStride, reinterpret_cast<void*>(3 * sizeof(float)));
    gl33::glEnableVertexAttribArray(1);
    gl33::glBindVertexArray(0);

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.08f, 0.08f, 0.1f, 1.0f);

    ShowWindow(hwnd, SW_SHOW);
    glViewport(0, 0, g_width, g_height);

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

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const double aspect = g_height > 0 ? static_cast<double>(g_width) / g_height : 1.0;
        const double nearPlane = std::max(0.01 * g_camera.distance, 1e-3);
        const double farPlane = g_camera.distance * 100.0 + 10.0;
        const Matrix4x4 projection = Matrix4x4::perspective(45.0 * 3.14159265358979323846 / 180.0, aspect,
                                                              nearPlane, farPlane);
        const Vector3 eye = g_camera.eye();
        const Matrix4x4 view = Matrix4x4::lookAt(eye, g_camera.center, Vector3(0.0, 1.0, 0.0));
        const Matrix4x4 mvp = projection * view;

        if (surfaceVertexCount > 0) {
            gl33::glUseProgram(surfaceProgram);
            gl33::glUniformMatrix4fv(surfaceMvpLoc, 1, GL_FALSE, mvp.data());
            gl33::glUniform3f(lightDirLoc, 1.0f, 1.0f, 1.0f);
            gl33::glUniform3f(baseColorLoc, 0.85f, 0.55f, 0.2f);
            gl33::glBindVertexArray(surfaceVao);
            glDrawArrays(GL_TRIANGLES, 0, surfaceVertexCount);
            gl33::glBindVertexArray(0);
        }

        gl33::glUseProgram(lineProgram);
        gl33::glUniformMatrix4fv(lineMvpLoc, 1, GL_FALSE, mvp.data());
        gl33::glBindVertexArray(lineVao);
        glDrawArrays(GL_LINES, 0, boxVertexCount);
        gl33::glBindVertexArray(0);

        SwapBuffers(hdc);
    }

    gl33::glDeleteBuffers(1, &surfaceVbo);
    gl33::glDeleteVertexArrays(1, &surfaceVao);
    gl33::glDeleteProgram(surfaceProgram);
    gl33::glDeleteBuffers(1, &lineVbo);
    gl33::glDeleteVertexArrays(1, &lineVao);
    gl33::glDeleteProgram(lineProgram);
    wglMakeCurrent(nullptr, nullptr);
    wglDeleteContext(glContext);
    ReleaseDC(hwnd, hdc);
    DestroyWindow(hwnd);
    return 0;
}

} // namespace isosurface_mode

// ===========================================================================
// Mode: turbulence (formerly apps/turbulence_viewer) -- u+ vs ln(y+).
// ===========================================================================
namespace turbulence_mode {

using aether::core::Matrix4x4;
using aether::solver::KEpsilonChannelFlowSolver1D;
using aether::solver::KOmegaSSTChannelFlowSolver1D;
using aether::solver::MixingLengthChannelFlowSolver1D;

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

int g_width = 900;
int g_height = 600;

const double kKarman = 0.41;
const double kLogLawB = 5.0;

double logYPlus(double yPlus) { return std::log(std::max(yPlus, 1e-6)); }

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

struct DrawRange {
    GLenum primitive;
    GLint first;
    GLsizei count;
};

int run() {
    const double height = 2.0;
    const double nu = 5e-5;
    const double source = 0.003;

    MixingLengthChannelFlowSolver1D mixingLength(300, height, nu, source);
    const std::size_t mixingLengthIterations = mixingLength.solve();
    const double uTauMixingLength = mixingLength.frictionVelocity();

    KEpsilonChannelFlowSolver1D kEpsilon(36, height, nu, source);
    const std::size_t kEpsilonIterations = kEpsilon.solve(20000, 1e-10);
    const double uTauKEpsilon = kEpsilon.frictionVelocity();

    KOmegaSSTChannelFlowSolver1D kOmegaSST(36, height, nu, source);
    const std::size_t kOmegaSSTIterations = kOmegaSST.solve(20000, 1e-10);
    const double uTauKOmegaSST = kOmegaSST.frictionVelocity();

    const double reTau = uTauMixingLength * (height / 2.0) / nu;
    std::printf("Aether Unified Viewer [turbulence] - u+ vs ln(y+), Re_tau ~ %.0f\n", reTau);
    std::printf("  comprimento de mistura: %zu iteracoes, u_tau=%.5f\n", mixingLengthIterations,
                uTauMixingLength);
    std::printf("  k-epsilon:              %zu iteracoes, u_tau=%.5f\n", kEpsilonIterations, uTauKEpsilon);
    std::printf("  k-omega SST:            %zu iteracoes, u_tau=%.5f\n", kOmegaSSTIterations,
                uTauKOmegaSST);

    const double xMin = logYPlus(1.0);
    const double xMax = logYPlus(2000.0);
    const double yMin = 0.0;
    const double yMax = 26.0;

    HWND hwnd = gl33::createSimpleWindow(L"AetherUnifiedViewerTurbulenceWindowClass",
                                          L"Aether Viewer - Turbulence", g_width, g_height, WndProc);
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

    std::vector<float> vertexData;
    std::vector<DrawRange> ranges;
    auto beginRange = [&]() { return static_cast<GLint>(vertexData.size() / 5); };
    auto endRange = [&](GLenum primitive, GLint first) {
        ranges.push_back({primitive, first, static_cast<GLsizei>(vertexData.size() / 5 - first)});
    };

    {
        const float c[3] = {0.6f, 0.6f, 0.65f};
        const GLint first = beginRange();
        pushVertex(vertexData, xMin, yMin, c[0], c[1], c[2]);
        pushVertex(vertexData, xMax, yMin, c[0], c[1], c[2]);
        pushVertex(vertexData, xMin, yMin, c[0], c[1], c[2]);
        pushVertex(vertexData, xMin, yMax, c[0], c[1], c[2]);
        endRange(GL_LINES, first);
    }
    {
        const float c[3] = {0.25f, 0.25f, 0.3f};
        const GLint first = beginRange();
        for (double decade : {1.0, 10.0, 100.0, 1000.0}) {
            const double x = logYPlus(decade);
            pushVertex(vertexData, x, yMin, c[0], c[1], c[2]);
            pushVertex(vertexData, x, yMax, c[0], c[1], c[2]);
        }
        for (double u = 0.0; u <= yMax; u += 5.0) {
            pushVertex(vertexData, xMin, u, c[0], c[1], c[2]);
            pushVertex(vertexData, xMax, u, c[0], c[1], c[2]);
        }
        endRange(GL_LINES, first);
    }
    {
        const float c[3] = {0.75f, 0.75f, 0.78f};
        const GLint first = beginRange();
        for (double yPlus = 1.0; yPlus <= 5.0; yPlus += 0.5) {
            pushVertex(vertexData, logYPlus(yPlus), yPlus, c[0], c[1], c[2]);
        }
        endRange(GL_LINE_STRIP, first);
    }
    {
        const float c[3] = {0.75f, 0.75f, 0.78f};
        const GLint first = beginRange();
        for (double yPlus = 30.0; yPlus <= 2000.0; yPlus *= 1.05) {
            const double uPlus = (1.0 / kKarman) * std::log(yPlus) + kLogLawB;
            pushVertex(vertexData, logYPlus(yPlus), uPlus, c[0], c[1], c[2]);
        }
        endRange(GL_LINE_STRIP, first);
    }
    {
        const float c[3] = {0.3f, 0.55f, 1.0f};
        const GLint first = beginRange();
        for (std::size_t j = 0; j < 150; ++j) {
            const double yPlus = mixingLength.wallDistance(j) * uTauMixingLength / nu;
            const double uPlus = mixingLength.u(j) / uTauMixingLength;
            pushVertex(vertexData, logYPlus(yPlus), uPlus, c[0], c[1], c[2]);
        }
        endRange(GL_LINE_STRIP, first);
    }
    {
        const float c[3] = {1.0f, 0.6f, 0.2f};
        const GLint firstLine = beginRange();
        for (std::size_t j = 0; j < 18; ++j) {
            const double yPlus = kEpsilon.wallDistance(j) * uTauKEpsilon / nu;
            const double uPlus = kEpsilon.u(j) / uTauKEpsilon;
            pushVertex(vertexData, logYPlus(yPlus), uPlus, c[0], c[1], c[2]);
        }
        endRange(GL_LINE_STRIP, firstLine);
        const GLint firstPoints = beginRange();
        for (std::size_t j = 0; j < 18; ++j) {
            const double yPlus = kEpsilon.wallDistance(j) * uTauKEpsilon / nu;
            const double uPlus = kEpsilon.u(j) / uTauKEpsilon;
            pushVertex(vertexData, logYPlus(yPlus), uPlus, c[0], c[1], c[2]);
        }
        endRange(GL_POINTS, firstPoints);
    }
    {
        const float c[3] = {0.3f, 0.9f, 0.4f};
        const GLint firstLine = beginRange();
        for (std::size_t j = 0; j < 18; ++j) {
            const double yPlus = kOmegaSST.wallDistance(j) * uTauKOmegaSST / nu;
            const double uPlus = kOmegaSST.u(j) / uTauKOmegaSST;
            pushVertex(vertexData, logYPlus(yPlus), uPlus, c[0], c[1], c[2]);
        }
        endRange(GL_LINE_STRIP, firstLine);
        const GLint firstPoints = beginRange();
        for (std::size_t j = 0; j < 18; ++j) {
            const double yPlus = kOmegaSST.wallDistance(j) * uTauKOmegaSST / nu;
            const double uPlus = kOmegaSST.u(j) / uTauKOmegaSST;
            pushVertex(vertexData, logYPlus(yPlus), uPlus, c[0], c[1], c[2]);
        }
        endRange(GL_POINTS, firstPoints);
    }

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

    glClearColor(0.08f, 0.08f, 0.1f, 1.0f);
    glLineWidth(2.0f);
    glPointSize(6.0f);

    ShowWindow(hwnd, SW_SHOW);
    glViewport(0, 0, g_width, g_height);

    const double marginFrac = 0.08;
    const double xSpan = xMax - xMin;
    const double ySpan = yMax - yMin;
    const Matrix4x4 projection = Matrix4x4::ortho(xMin - marginFrac * xSpan, xMax + marginFrac * xSpan,
                                                   yMin - marginFrac * ySpan, yMax + marginFrac * ySpan,
                                                   -1.0, 1.0);

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
        for (const DrawRange& range : ranges) {
            glDrawArrays(range.primitive, range.first, range.count);
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

} // namespace turbulence_mode

// ===========================================================================
// Entry point: mode dispatch.
// ===========================================================================

// ===========================================================================
// Mode: sim (Module 9) -- live 2D lid-driven cavity driven from a CAD-like
// control panel: pick the turbulence closure, edit parameters, run/pause/
// step/reset, and read diagnostics, all without recompiling.
//
// **This is the first mode that is interactive in the CFD sense** rather
// than the camera sense: every earlier mode solved its problem up front and
// then displayed a fixed result. Here the solver advances inside the render
// loop, and the panel writes back into it.
//
// Two deliberate simplifications, both to keep this first pass about the UI
// rather than about renderer plumbing:
//   - The field is drawn as a per-cell speed heatmap through the UI layer's
//     own drawRect(), not a separate GL pipeline. The UI already batches
//     screen-space coloured rectangles into a single draw call, so a 40x40
//     cavity is ~1600 more rectangles in the same batch -- cheaper and far
//     less code than standing up a second shader, and it exercises the new
//     UI layer hard enough to double as its smoke test.
//   - 2D, not 3D. The panel and its wiring are dimension-agnostic, so
//     pointing them at the 3D staggered cavity is mechanical follow-up --
//     but the 3D closures cost enough per step that interactive stepping
//     would need a worker thread, which is its own task.
// ===========================================================================
namespace sim_mode {

HWND g_window = nullptr;
int g_width = 1280;
int g_height = 800;
gl33::UiInput g_input;
bool g_pendingPress = false;
bool g_pendingRelease = false;
// WM_CHAR characters queued since the last frame collected them into
// g_input.textInput -- same "queue in WndProc, drain once per frame"
// pattern as g_pendingPress/g_pendingRelease above, needed because WndProc
// can fire multiple times between two iterations of the render loop.
std::string g_pendingText;

LRESULT CALLBACK WndProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam) {
    switch (message) {
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    case WM_SIZE:
        g_width = LOWORD(lParam);
        g_height = HIWORD(lParam);
        glViewport(0, 0, g_width, g_height);
        return 0;
    case WM_MOUSEMOVE:
        g_input.mouseX = GET_X_LPARAM(lParam);
        g_input.mouseY = GET_Y_LPARAM(lParam);
        return 0;
    case WM_LBUTTONDOWN:
        g_input.mouseDown = true;
        g_pendingPress = true;
        SetCapture(hwnd);
        return 0;
    case WM_LBUTTONUP:
        g_input.mouseDown = false;
        g_pendingRelease = true;
        ReleaseCapture();
        return 0;
    case WM_CHAR:
        // TranslateMessage() (already called in this app's message loop)
        // turns WM_KEYDOWN into WM_CHAR, including for Backspace/Enter/Esc
        // -- see UiInput::textInput's own comment for why that is enough.
        if (wParam > 0 && wParam < 0x80) {
            g_pendingText.push_back(static_cast<char>(wParam));
        }
        return 0;
    default:
        return DefWindowProc(hwnd, message, wParam, lParam);
    }
}

// The four 2D cavity closures share an identical constructor signature and
// an identical u/v/time/maxDivergence surface, so holding one owning pointer
// per closure plus a small visitor is enough -- no wrapper hierarchy needed.
enum class Closure { Laminar = 0, MixingLength, KEpsilon, KOmegaSST };

const char* closureName(Closure closure) {
    switch (closure) {
    case Closure::Laminar:
        return "laminar";
    case Closure::MixingLength:
        return "mixing-length";
    case Closure::KEpsilon:
        return "k-epsilon";
    case Closure::KOmegaSST:
        return "k-omega SST";
    }
    return "?";
}

struct Simulation {
    Closure closure = Closure::Laminar;
    std::size_t n = 40;
    double viscosity = 0.01;
    double lidVelocity = 1.0;
    long long steps = 0;

    std::unique_ptr<aether::solver::LidDrivenCavitySolver2D> laminar;
    std::unique_ptr<aether::solver::MixingLengthLidDrivenCavitySolver2D> mixingLength;
    std::unique_ptr<aether::solver::KEpsilonLidDrivenCavitySolver2D> kEpsilon;
    std::unique_ptr<aether::solver::KOmegaSSTLidDrivenCavitySolver2D> kOmegaSST;

    void rebuild() {
        laminar.reset();
        mixingLength.reset();
        kEpsilon.reset();
        kOmegaSST.reset();
        steps = 0;
        switch (closure) {
        case Closure::Laminar:
            laminar = std::make_unique<aether::solver::LidDrivenCavitySolver2D>(n, n, 1.0, 1.0, viscosity,
                                                                                 lidVelocity);
            break;
        case Closure::MixingLength:
            mixingLength = std::make_unique<aether::solver::MixingLengthLidDrivenCavitySolver2D>(
                n, n, 1.0, 1.0, viscosity, lidVelocity);
            break;
        case Closure::KEpsilon:
            kEpsilon = std::make_unique<aether::solver::KEpsilonLidDrivenCavitySolver2D>(
                n, n, 1.0, 1.0, viscosity, lidVelocity);
            break;
        case Closure::KOmegaSST:
            kOmegaSST = std::make_unique<aether::solver::KOmegaSSTLidDrivenCavitySolver2D>(
                n, n, 1.0, 1.0, viscosity, lidVelocity);
            break;
        }
    }

    template <typename Fn> auto visit(Fn&& fn) const -> decltype(fn(*laminar)) {
        if (laminar) {
            return fn(*laminar);
        }
        if (mixingLength) {
            return fn(*mixingLength);
        }
        if (kEpsilon) {
            return fn(*kEpsilon);
        }
        return fn(*kOmegaSST);
    }

    void step() {
        const auto advance = [](auto& solver) { solver.step(0.3 * solver.stableTimeStep()); };
        if (laminar) {
            advance(*laminar);
        } else if (mixingLength) {
            advance(*mixingLength);
        } else if (kEpsilon) {
            advance(*kEpsilon);
        } else {
            advance(*kOmegaSST);
        }
        ++steps;
    }

    double speedAt(std::size_t i, std::size_t j) const {
        return visit([&](const auto& solver) {
            const double u = solver.u(i, j);
            const double v = solver.v(i, j);
            return std::sqrt(u * u + v * v);
        });
    }
    double pressureAt(std::size_t i, std::size_t j) const {
        return visit([&](const auto& solver) { return solver.pressure(i, j); });
    }
    double maxDivergence() const {
        return visit([](const auto& solver) { return solver.maxDivergence(); });
    }
    double time() const {
        return visit([](const auto& solver) { return solver.time(); });
    }
};

int run() {
    std::printf("Aether Unified Viewer - modo 'sim' (Modulo 9: UI com paineis)\n");
    std::printf("  painel a esquerda: fechamento, parametros, run/pause/passo/reiniciar\n");
    std::printf("  o campo e' redesenhado a cada quadro a partir do solver ao vivo\n\n");
    std::fflush(stdout);

    g_window = gl33::createSimpleWindow(L"AetherSim", L"Aether - Simulacao interativa (Modulo 9)", g_width,
                                        g_height, WndProc);
    if (g_window == nullptr) {
        std::fprintf(stderr, "erro ao criar janela\n");
        return 1;
    }
    HDC hdc = nullptr;
    HGLRC context = gl33::createGl33Context(g_window, hdc);
    if (context == nullptr) {
        return 1;
    }

    gl33::Ui ui;
    if (!ui.initialize()) {
        std::fprintf(stderr, "erro ao inicializar a UI\n");
        return 1;
    }

    Simulation sim;
    sim.rebuild();

    // Starts running: opening this mode to a frozen, uniformly-blue field
    // reads as broken, and the first thing anyone would do is press RODAR.
    bool running = true;
    bool showPressure = false;
    double resolution = static_cast<double>(sim.n);
    double viscosity = sim.viscosity;
    double lidVelocity = sim.lidVelocity;
    constexpr int kStepsPerFrame = 4;

    ShowWindow(g_window, SW_SHOW);
    UpdateWindow(g_window);

    MSG message{};
    bool quit = false;
    while (!quit) {
        while (PeekMessage(&message, nullptr, 0, 0, PM_REMOVE)) {
            if (message.message == WM_QUIT) {
                quit = true;
            }
            TranslateMessage(&message);
            DispatchMessage(&message);
        }
        if (quit) {
            break;
        }

        g_input.mousePressed = g_pendingPress;
        g_input.mouseReleased = g_pendingRelease;
        g_input.textInput = g_pendingText;
        g_pendingPress = false;
        g_pendingRelease = false;
        g_pendingText.clear();

        if (running) {
            for (int s = 0; s < kStepsPerFrame; ++s) {
                sim.step();
            }
        }

        glClearColor(0.06f, 0.07f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        ui.begin(g_width, g_height, g_input);

        // --- field: per-cell heatmap, emitted into the UI's own batch ---
        const int panelWidth = 260;
        const int margin = 24;
        const int available = std::min(g_width - panelWidth - 3 * margin, g_height - 2 * margin);
        const int fieldSize = std::max(available, 64);
        const int fieldX = panelWidth + 2 * margin;
        const int fieldY = (g_height - fieldSize) / 2;

        double maxValue = 1e-12;
        double minValue = 0.0;
        for (std::size_t j = 0; j < sim.n; ++j) {
            for (std::size_t i = 0; i < sim.n; ++i) {
                const double value = showPressure ? sim.pressureAt(i, j) : sim.speedAt(i, j);
                maxValue = std::max(maxValue, value);
                minValue = std::min(minValue, value);
            }
        }
        const double span = (maxValue - minValue) > 1e-12 ? (maxValue - minValue) : 1.0;

        const int cell = std::max(fieldSize / static_cast<int>(sim.n), 1);
        for (std::size_t j = 0; j < sim.n; ++j) {
            for (std::size_t i = 0; i < sim.n; ++i) {
                const double value = showPressure ? sim.pressureAt(i, j) : sim.speedAt(i, j);
                float r = 0.0f;
                float g = 0.0f;
                float b = 0.0f;
                colorForValue((value - minValue) / span, r, g, b);
                // Row j=0 is the cavity floor, so flip j into screen space.
                const int px = fieldX + static_cast<int>(i) * cell;
                const int py = fieldY + (static_cast<int>(sim.n) - 1 - static_cast<int>(j)) * cell;
                ui.drawRect(px, py, cell, cell, gl33::UiColor{r, g, b, 1.0f});
            }
        }
        ui.drawText(fieldX, fieldY - 16,
                    showPressure ? "pressao (azul=baixa, vermelho=alta)"
                                 : "velocidade |u| (azul=baixa, vermelho=alta)",
                    gl33::UiColor{0.60f, 0.63f, 0.67f, 1.0f});

        // --- control panel ---
        ui.beginPanel(margin, margin, panelWidth, "AETHER - CAVIDADE 2D");

        ui.label("Fechamento de turbulencia");
        for (int c = 0; c < 4; ++c) {
            const auto candidate = static_cast<Closure>(c);
            std::string caption = (sim.closure == candidate) ? "[x] " : "[ ] ";
            caption += closureName(candidate);
            if (ui.button(caption) && sim.closure != candidate) {
                sim.closure = candidate;
                sim.rebuild();
                running = false;
            }
        }

        ui.separator();
        ui.label("Parametros");
        bool needsRebuild = false;
        // Resolution stays a slider: it must land on an integer and a
        // coarse drag is the natural way to pick one. Viscosity and lid
        // velocity get textField() instead -- both often need a precise
        // value (a specific Reynolds number) that a 260px slider track
        // cannot address finely, which is the whole reason M9.4 added
        // keyboard entry on top of the slider widget.
        needsRebuild |= ui.slider("resolucao", &resolution, 16, 80);
        needsRebuild |= ui.textField("viscosidade", &viscosity, 0.0005, 0.5);
        needsRebuild |= ui.textField("veloc. da tampa", &lidVelocity, 0.0, 5.0);
        if (needsRebuild) {
            sim.n = static_cast<std::size_t>(resolution);
            sim.viscosity = viscosity;
            sim.lidVelocity = lidVelocity;
            sim.rebuild();
            running = false;
        }

        ui.separator();
        ui.label("Simulacao");
        if (ui.button(running ? "PAUSAR" : "RODAR")) {
            running = !running;
        }
        if (ui.button("+1 PASSO")) {
            sim.step();
        }
        if (ui.button("REINICIAR")) {
            sim.rebuild();
            running = false;
        }
        ui.checkbox("mostrar pressao", &showPressure);

        ui.separator();
        ui.label("Diagnostico");
        char buffer[128];
        std::snprintf(buffer, sizeof(buffer), "Re = %.0f", sim.lidVelocity / sim.viscosity);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "malha: %zux%zu", sim.n, sim.n);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "passos: %lld", sim.steps);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "t = %.4f s", sim.time());
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "div max = %.2e", sim.maxDivergence());
        ui.label(buffer);

        ui.endPanel();
        ui.end();

        SwapBuffers(hdc);
        Sleep(1);
    }

    ui.shutdown();
    wglMakeCurrent(nullptr, nullptr);
    wglDeleteContext(context);
    ReleaseDC(g_window, hdc);
    return 0;
}

} // namespace sim_mode

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "uso: aether_unified_viewer <modo> [args]\n");
        std::fprintf(stderr, "  modos disponiveis:\n");
        std::fprintf(stderr, "    mesh <arquivo.stl>   visualizador de malha STL (Modulo 8)\n");
        std::fprintf(stderr, "    heatmap              conducao de calor 2D (Modulo 5/7)\n");
        std::fprintf(stderr, "    cavity                cavidade com tampa deslizante + streamlines (Modulo 4/7)\n");
        std::fprintf(stderr, "    cavity3d              cavidade 3D staggered, campo de vetores (Modulo 4/8)\n");
        std::fprintf(stderr, "    isosurface            iso-superficie de nu_t, k-omega SST 3D (Modulo 6/7/8)\n");
        std::fprintf(stderr, "    turbulence            u+ vs ln(y+) dos fechamentos de turbulencia (Modulo 6/7)\n");
        return 1;
    }

    const char* mode = argv[1];
    if (std::strcmp(mode, "mesh") == 0) {
        return mesh_mode::run(argc, argv);
    }
    if (std::strcmp(mode, "heatmap") == 0) {
        return heatmap_mode::run();
    }
    if (std::strcmp(mode, "cavity") == 0) {
        return cavity_mode::run();
    }
    if (std::strcmp(mode, "cavity3d") == 0) {
        return cavity3d_mode::run();
    }
    if (std::strcmp(mode, "isosurface") == 0) {
        return isosurface_mode::run();
    }
    if (std::strcmp(mode, "turbulence") == 0) {
        return turbulence_mode::run();
    }
    if (std::strcmp(mode, "sim") == 0) {
        return sim_mode::run();
    }

    std::fprintf(stderr, "modo desconhecido: %s\n", mode);
    std::fprintf(stderr, "modos disponiveis: mesh, heatmap, cavity, cavity3d, isosurface, turbulence\n");
    return 1;
}
