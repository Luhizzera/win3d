// Aether Viewer - Module 8: a Win32 + modern (OpenGL 3.3 core profile,
// shader/VBO-based) window that loads an STL file and displays it with an
// orbit camera. Uses the shared apps/common/Gl33 toolkit (context
// bootstrap, shader compile/link) rather than repeating that boilerplate
// here -- see Gl33.hpp for the full rationale (no GLFW/GLAD/Vulkan,
// hand-loaded GL 3.3 function pointers, the throwaway-context trick).
//
// Scope note: this is the first of the roadmap's four legacy viewer apps
// (viewer, field_viewer, cavity_viewer, turbulence_viewer) to move to the
// modern pipeline, deliberately -- it is the most complex 3D case (lit
// triangle mesh vs. the others' flat 2D overlays), so proving the shader/
// VBO/WGL-core-context machinery here first is lower-risk than attempting
// all four (and unifying them into one app) at once. The other three, and
// full unification into a single app, remain explicitly open follow-up
// work.

#include "aether/core/Matrix4x4.hpp"
#include "aether/geometry/StlIO.hpp"
#include "aether/geometry/TriangleMesh.hpp"
#include "aether_app/Gl33.hpp"

#include <windowsx.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using aether::core::Matrix4x4;
using aether::core::Vector3;
using aether::geometry::loadStl;
using aether::geometry::TriangleMesh;
namespace gl33 = aether::app;

namespace {

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

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "uso: aether_viewer <arquivo.stl>\n");
        return 1;
    }

    TriangleMesh mesh;
    try {
        mesh = loadStl(argv[1]);
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

    std::printf("Aether Viewer - %s\n", argv[1]);
    std::printf("  vertices: %zu, triangulos: %zu\n", mesh.vertexCount(), mesh.triangleCount());
    std::printf("  area: %.6f, volume: %.6f, watertight: %s\n", mesh.surfaceArea(), mesh.volume(),
                mesh.isWatertight() ? "sim" : "nao");

    HWND hwnd = gl33::createSimpleWindow(L"AetherViewerWindowClass", L"Aether Viewer", g_width, g_height,
                                          WndProc);
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

    // Flat per-face shading: each triangle contributes 3 vertices with its
    // own face normal repeated, matching the old glNormal-before-glVertex
    // immediate-mode behavior exactly (no shared/smoothed vertex normals).
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
        const Matrix4x4 mvp = projection * view; // model is identity

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
