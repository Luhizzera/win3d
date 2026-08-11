// Aether Turbulence Viewer - Module 6/7: the classic u+ vs ln(y+) semi-log
// plot used throughout turbulence modeling to check a closure against the
// log law. Solves the same fully-developed turbulent channel flow with
// MixingLengthChannelFlowSolver1D (fine resolution, resolved to the wall),
// KEpsilonChannelFlowSolver1D, and KOmegaSSTChannelFlowSolver1D (both
// coarser wall-function meshes), and overlays the theoretical log law
// (log-law region) and viscous sublayer (u+ = y+) lines solver_tests.cpp
// validates all three models against. Modernized (fourth of the four
// legacy viewer apps) onto the shared OpenGL 3.3 core-profile / shader /
// VBO pipeline in apps/common/Gl33.
//
// No text rendering (no font system in this minimal GL setup, same as the
// other viewer apps) -- axis ranges and a legend are printed to the
// console instead. Gridlines at round y+ decades (1, 10, 100, 1000) and
// u+ values stand in for axis labels.

#include "aether/core/Matrix4x4.hpp"
#include "aether/solver/KEpsilonChannelFlowSolver1D.hpp"
#include "aether/solver/KOmegaSSTChannelFlowSolver1D.hpp"
#include "aether/solver/MixingLengthChannelFlowSolver1D.hpp"
#include "aether_app/Gl33.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using aether::core::Matrix4x4;
using aether::solver::KEpsilonChannelFlowSolver1D;
using aether::solver::KOmegaSSTChannelFlowSolver1D;
using aether::solver::MixingLengthChannelFlowSolver1D;
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

void pushVertex(std::vector<float>& out, double x, double y, float r, float g, float b) {
    out.push_back(static_cast<float>(x));
    out.push_back(static_cast<float>(y));
    out.push_back(r);
    out.push_back(g);
    out.push_back(b);
}

// One contiguous run of vertices in the shared VBO, drawn with a single
// glDrawArrays call.
struct DrawRange {
    GLenum primitive;
    GLint first;
    GLsizei count;
};

} // namespace

int main() {
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
    std::printf("Aether Turbulence Viewer - u+ vs ln(y+), Re_tau ~ %.0f\n", reTau);
    std::printf("  comprimento de mistura: %zu iteracoes, u_tau=%.5f\n", mixingLengthIterations,
                uTauMixingLength);
    std::printf("  k-epsilon:              %zu iteracoes, u_tau=%.5f\n", kEpsilonIterations, uTauKEpsilon);
    std::printf("  k-omega SST:            %zu iteracoes, u_tau=%.5f\n", kOmegaSSTIterations,
                uTauKOmegaSST);
    std::printf("  linhas: cinza tracejada = lei log teorica (1/kappa * ln(y+) + %.1f) e u+=y+\n",
                kLogLawB);
    std::printf("          azul = comprimento de mistura (resolvido ate a parede)\n");
    std::printf("          laranja = k-epsilon (funcoes de parede, y+ >= ~30)\n");
    std::printf("          verde = k-omega SST (funcoes de parede, y+ >= ~30)\n\n");
    std::printf("  esperado: a curva laranja (k-epsilon) fica colada na reta teorica desde o\n");
    std::printf("  primeiro ponto - a celula adjacente a parede e FIXADA por essa mesma formula\n");
    std::printf("  (e' a condicao de contorno, nao uma coincidencia). Ja a curva azul (comprimento\n");
    std::printf("  de mistura) fica sistematicamente ABAIXO da reta, com a MESMA INCLINACAO -\n");
    std::printf("  o modelo acerta o expoente 1/kappa (validado em solver_tests.cpp) mas nao foi\n");
    std::printf("  ajustado pra reproduzir a constante aditiva B=5.0 (essa constante nao existe\n");
    std::printf("  explicitamente nesse modelo mais simples, sem amortecimento de van Driest);\n");
    std::printf("  isso e' esperado, nao um bug.\n\n");

    const double xMin = logYPlus(1.0);
    const double xMax = logYPlus(2000.0);
    const double yMin = 0.0;
    const double yMax = 26.0;

    HWND hwnd = gl33::createSimpleWindow(L"AetherTurbulenceViewerWindowClass", L"Aether Turbulence Viewer",
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

    std::vector<float> vertexData;
    std::vector<DrawRange> ranges;
    auto beginRange = [&]() { return static_cast<GLint>(vertexData.size() / 5); };
    auto endRange = [&](GLenum primitive, GLint first) {
        ranges.push_back({primitive, first, static_cast<GLsizei>(vertexData.size() / 5 - first)});
    };

    // Axes.
    {
        const float c[3] = {0.6f, 0.6f, 0.65f};
        const GLint first = beginRange();
        pushVertex(vertexData, xMin, yMin, c[0], c[1], c[2]);
        pushVertex(vertexData, xMax, yMin, c[0], c[1], c[2]);
        pushVertex(vertexData, xMin, yMin, c[0], c[1], c[2]);
        pushVertex(vertexData, xMin, yMax, c[0], c[1], c[2]);
        endRange(GL_LINES, first);
    }
    // Gridlines.
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
    // Theoretical viscous sublayer (u+ = y+, y+ in [1,5]).
    {
        const float c[3] = {0.75f, 0.75f, 0.78f};
        const GLint first = beginRange();
        for (double yPlus = 1.0; yPlus <= 5.0; yPlus += 0.5) {
            pushVertex(vertexData, logYPlus(yPlus), yPlus, c[0], c[1], c[2]);
        }
        endRange(GL_LINE_STRIP, first);
    }
    // Theoretical log law (y+ in [30, 2000]).
    {
        const float c[3] = {0.75f, 0.75f, 0.78f};
        const GLint first = beginRange();
        for (double yPlus = 30.0; yPlus <= 2000.0; yPlus *= 1.05) {
            const double uPlus = (1.0 / kKarman) * std::log(yPlus) + kLogLawB;
            pushVertex(vertexData, logYPlus(yPlus), uPlus, c[0], c[1], c[2]);
        }
        endRange(GL_LINE_STRIP, first);
    }
    // Mixing length: resolved from wall to centerline.
    {
        const float c[3] = {0.3f, 0.55f, 1.0f};
        const GLint first = beginRange();
        for (std::size_t j = 0; j < 150; ++j) { // first half; symmetric about centerline
            const double yPlus = mixingLength.wallDistance(j) * uTauMixingLength / nu;
            const double uPlus = mixingLength.u(j) / uTauMixingLength;
            pushVertex(vertexData, logYPlus(yPlus), uPlus, c[0], c[1], c[2]);
        }
        endRange(GL_LINE_STRIP, first);
    }
    // k-epsilon: coarser wall-function mesh (line + point markers).
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
    // k-omega SST: same coarser wall-function mesh as k-epsilon.
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
