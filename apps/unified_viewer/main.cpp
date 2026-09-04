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
#include "aether/solver/DesSstLidDrivenCavitySolver3D.hpp"
#include "aether/solver/KEpsilonLidDrivenCavitySolver3D.hpp"
#include "aether/solver/KOmegaSSTLidDrivenCavitySolver3D.hpp"
#include "aether/solver/MixingLengthLidDrivenCavitySolver3D.hpp"
#include "aether/solver/SmagorinskyLesLidDrivenCavitySolver3D.hpp"
#include "aether/solver/StaggeredLidDrivenCavitySolver3D.hpp"
#include "aether/solver/SteadyDiffusionSolver.hpp"
#include "aether_app/Gl33.hpp"
#include "aether_app/Ui.hpp"
#include "Mode.hpp"

namespace aether_app_modes {
// Shared layout constants: the Workspace's own mode-switcher panel always
// occupies the left column (see Workspace below), so any mode that draws
// its own ui.beginPanel() -- currently SimMode and Sim3DMode -- must start
// past it rather than at the same x=kSidebarMargin the switcher itself
// uses, or the two panels would draw on top of each other (Ui supports
// only one open panel at a time, with no auto-flow between separate
// beginPanel() calls -- see Ui.hpp's own comment on that).
constexpr int kSidebarMargin = 24;
constexpr int kModeSwitcherWidth = 220;
constexpr int kModePanelX = kSidebarMargin * 2 + kModeSwitcherWidth;
} // namespace aether_app_modes

#include <windowsx.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
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

namespace aether_app_modes {

class MeshMode : public Mode {
public:
    explicit MeshMode(std::string stlPath = {}) : stlPath_(std::move(stlPath)) {}

    const char* label() const override { return "MALHA (mesh)"; }

    bool init() override {
        using aether::core::Matrix4x4;
        using aether::core::Vector3;
        using aether::geometry::loadStl;
        using aether::geometry::TriangleMesh;

        TriangleMesh mesh;
        if (!stlPath_.empty()) {
            try {
                mesh = loadStl(stlPath_);
            } catch (const std::exception& e) {
                std::fprintf(stderr, "erro ao carregar STL: %s\n", e.what());
                return false;
            }
            std::printf("Aether Workspace [mesh] - %s\n", stlPath_.c_str());
        } else {
            mesh = buildIcosahedron();
            std::printf("Aether Workspace [mesh] - icosaedro construido em memoria "
                        "(passe um .stl na linha de comando para trocar)\n");
        }

        std::printf("  vertices: %zu, triangulos: %zu\n", mesh.vertexCount(), mesh.triangleCount());
        std::printf("  area: %.6f, volume: %.6f, watertight: %s\n", mesh.surfaceArea(), mesh.volume(),
                    mesh.isWatertight() ? "sim" : "nao");

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
        camera_.center = (minBound + maxBound) * 0.5;
        const Vector3 extent = maxBound - minBound;
        const double diagonal = std::sqrt(extent.x * extent.x + extent.y * extent.y + extent.z * extent.z);
        camera_.distance = diagonal > 0 ? diagonal * 1.5 : 5.0;

        const GLuint vertexShader = gl33::compileShader(gl33::kGlVertexShader, kVertexShaderSource);
        const GLuint fragmentShader = gl33::compileShader(gl33::kGlFragmentShader, kFragmentShaderSource);
        program_ = gl33::linkProgram(vertexShader, fragmentShader);
        gl33::glDeleteShader(vertexShader);
        gl33::glDeleteShader(fragmentShader);
        mvpLoc_ = gl33::glGetUniformLocation(program_, "uMvp");
        lightDirLoc_ = gl33::glGetUniformLocation(program_, "uLightDir");
        baseColorLoc_ = gl33::glGetUniformLocation(program_, "uBaseColor");

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
        vertexCount_ = static_cast<GLsizei>(mesh.triangleCount() * 3);

        gl33::glGenVertexArrays(1, &vao_);
        gl33::glBindVertexArray(vao_);
        gl33::glGenBuffers(1, &vbo_);
        gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo_);
        gl33::glBufferData(gl33::kGlArrayBuffer,
                            static_cast<gl33::GLsizeiptr>(vertexData.size() * sizeof(float)), vertexData.data(),
                            gl33::kGlStaticDraw);
        const GLsizei stride = 6 * sizeof(float);
        gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
        gl33::glEnableVertexAttribArray(0);
        gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));
        gl33::glEnableVertexAttribArray(1);
        gl33::glBindVertexArray(0);
        return true;
    }

    void shutdown() override {
        gl33::glDeleteBuffers(1, &vbo_);
        gl33::glDeleteVertexArrays(1, &vao_);
        gl33::glDeleteProgram(program_);
        dragging_ = false;
    }

    void handleInput(const aether::app::UiInput& input, bool uiWantsMouse, double wheelDelta) override {
        if (wheelDelta != 0.0) {
            const double factor = wheelDelta > 0 ? 0.9 : 1.1;
            camera_.distance = std::max(1e-3, camera_.distance * factor);
        }
        if (input.mousePressed && !uiWantsMouse) {
            dragging_ = true;
            lastMouseX_ = input.mouseX;
            lastMouseY_ = input.mouseY;
        }
        if (input.mouseReleased) {
            dragging_ = false;
        }
        if (dragging_ && input.mouseDown) {
            camera_.yaw += (input.mouseX - lastMouseX_) * 0.01;
            camera_.pitch = std::clamp(camera_.pitch + (input.mouseY - lastMouseY_) * 0.01, -1.5, 1.5);
            lastMouseX_ = input.mouseX;
            lastMouseY_ = input.mouseY;
        }
    }

    void update() override {}

    void renderScene(int width, int height) override {
        using aether::core::Matrix4x4;
        using aether::core::Vector3;

        glEnable(GL_DEPTH_TEST);
        glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const double aspect = height > 0 ? static_cast<double>(width) / height : 1.0;
        const double nearPlane = std::max(0.01 * camera_.distance, 1e-3);
        const double farPlane = camera_.distance * 100.0 + 10.0;
        const Matrix4x4 projection = Matrix4x4::perspective(45.0 * 3.14159265358979323846 / 180.0, aspect,
                                                              nearPlane, farPlane);
        const Vector3 eye = camera_.eye();
        const Matrix4x4 view = Matrix4x4::lookAt(eye, camera_.center, Vector3(0.0, 1.0, 0.0));
        const Matrix4x4 mvp = projection * view;

        gl33::glUseProgram(program_);
        gl33::glUniformMatrix4fv(mvpLoc_, 1, GL_FALSE, mvp.data());
        gl33::glUniform3f(lightDirLoc_, 1.0f, 1.0f, 1.0f);
        gl33::glUniform3f(baseColorLoc_, 0.7f, 0.75f, 0.85f);

        gl33::glBindVertexArray(vao_);
        glDrawArrays(GL_TRIANGLES, 0, vertexCount_);
        gl33::glBindVertexArray(0);
    }

    void renderPanel(aether::app::Ui&) override {}

private:
    struct OrbitCamera {
        aether::core::Vector3 center;
        double distance = 5.0;
        double yaw = 0.7;
        double pitch = 0.5;

        aether::core::Vector3 eye() const {
            const double cp = std::cos(pitch);
            return center + aether::core::Vector3(distance * cp * std::cos(yaw), distance * std::sin(pitch),
                                                    distance * cp * std::sin(yaw));
        }
    };

    static aether::geometry::TriangleMesh buildIcosahedron() {
        const double phi = (1.0 + std::sqrt(5.0)) / 2.0;
        const double raw[12][3] = {{-1, phi, 0}, {1, phi, 0},  {-1, -phi, 0}, {1, -phi, 0},
                                   {0, -1, phi}, {0, 1, phi},  {0, -1, -phi}, {0, 1, -phi},
                                   {phi, 0, -1}, {phi, 0, 1},  {-phi, 0, -1}, {-phi, 0, 1}};
        const std::size_t faces[20][3] = {
            {0, 11, 5}, {0, 5, 1},  {0, 1, 7},  {0, 7, 10}, {0, 10, 11}, {1, 5, 9},  {5, 11, 4}, {11, 10, 2},
            {10, 7, 6}, {7, 1, 8},  {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},  {4, 9, 5},
            {2, 4, 11}, {6, 2, 10}, {8, 6, 7},  {9, 8, 1}};
        const double norm = std::sqrt(1.0 + phi * phi);

        aether::geometry::TriangleMesh mesh;
        for (const auto& v : raw) {
            mesh.addVertex(aether::core::Vector3(v[0] / norm, v[1] / norm, v[2] / norm));
        }
        for (const auto& f : faces) {
            mesh.addTriangle(f[0], f[1], f[2]);
        }
        return mesh;
    }

    static constexpr const char* kVertexShaderSource = R"GLSL(
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

    static constexpr const char* kFragmentShaderSource = R"GLSL(
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

    std::string stlPath_;
    OrbitCamera camera_;
    bool dragging_ = false;
    int lastMouseX_ = 0;
    int lastMouseY_ = 0;

    GLuint program_ = 0;
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLint mvpLoc_ = -1;
    GLint lightDirLoc_ = -1;
    GLint baseColorLoc_ = -1;
    GLsizei vertexCount_ = 0;
};

} // namespace aether_app_modes

// ===========================================================================
// Mode: heatmap (formerly apps/field_viewer) -- steady heat conduction.
// ===========================================================================

namespace aether_app_modes {

class HeatmapMode : public Mode {
public:
    const char* label() const override { return "CONDUCAO DE CALOR (heatmap)"; }

    bool init() override {
        using aether::core::Matrix4x4;
        using aether::core::Vector3;
        using aether::mesh::StructuredGrid3D;
        using aether::solver::SteadyDiffusionSolver;

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

        std::printf("Aether Workspace [heatmap] - conducao de calor estacionaria em grade %zux%zu\n", nx, ny);
        std::printf("  iteracoes ate convergir: %zu\n", iterations);

        const GLuint vertexShader = gl33::compileShader(gl33::kGlVertexShader, kVertexShaderSource);
        const GLuint fragmentShader = gl33::compileShader(gl33::kGlFragmentShader, kFragmentShaderSource);
        program_ = gl33::linkProgram(vertexShader, fragmentShader);
        gl33::glDeleteShader(vertexShader);
        gl33::glDeleteShader(fragmentShader);
        projectionLoc_ = gl33::glGetUniformLocation(program_, "uProjection");

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
        vertexCount_ = static_cast<GLsizei>(nx * ny * 6);

        gl33::glGenVertexArrays(1, &vao_);
        gl33::glBindVertexArray(vao_);
        gl33::glGenBuffers(1, &vbo_);
        gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo_);
        gl33::glBufferData(gl33::kGlArrayBuffer,
                            static_cast<gl33::GLsizeiptr>(vertexData.size() * sizeof(float)), vertexData.data(),
                            gl33::kGlStaticDraw);
        const GLsizei stride = 5 * sizeof(float);
        gl33::glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
        gl33::glEnableVertexAttribArray(0);
        gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(2 * sizeof(float)));
        gl33::glEnableVertexAttribArray(1);
        gl33::glBindVertexArray(0);

        projection_ = Matrix4x4::ortho(0.0, lx, 0.0, ly, -1.0, 1.0);
        return true;
    }

    void shutdown() override {
        gl33::glDeleteBuffers(1, &vbo_);
        gl33::glDeleteVertexArrays(1, &vao_);
        gl33::glDeleteProgram(program_);
    }

    void handleInput(const aether::app::UiInput&, bool, double) override {}
    void update() override {}

    void renderScene(int, int) override {
        glClearColor(0.05f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        gl33::glUseProgram(program_);
        gl33::glUniformMatrix4fv(projectionLoc_, 1, GL_FALSE, projection_.data());
        gl33::glBindVertexArray(vao_);
        glDrawArrays(GL_TRIANGLES, 0, vertexCount_);
        gl33::glBindVertexArray(0);
    }

    void renderPanel(aether::app::Ui&) override {}

private:
    static constexpr const char* kVertexShaderSource = R"GLSL(
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

    static constexpr const char* kFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
}
)GLSL";

    GLuint program_ = 0;
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLint projectionLoc_ = -1;
    GLsizei vertexCount_ = 0;
    aether::core::Matrix4x4 projection_;
};

} // namespace aether_app_modes

// ===========================================================================
// Mode: cavity (formerly apps/cavity_viewer) -- velocity arrows + streamlines.
// ===========================================================================

// ===========================================================================
// Workspace-hosted version of the mode above. See Mode.hpp for the
// interface contract and apps/unified_viewer's migration plan for why this
// exists alongside (temporarily) the original namespace: the original
// stays until every one of the 8 modes has an equivalent here and the new
// Workspace (added once all 8 exist) is verified end to end, at which
// point every `*_mode` namespace above is deleted in one pass.
namespace aether_app_modes {

class CavityMode : public Mode {
public:
    const char* label() const override { return "CAVIDADE 2D (cavity)"; }

    bool init() override {
        using aether::core::Matrix4x4;
        using aether::solver::LidDrivenCavitySolver2D;
        using aether::postprocessing::Streamline2D;

        const std::size_t n = 24;
        const double length = 1.0;
        const double viscosity = 0.1;
        const double lidVelocity = 1.0; // Re = lidVelocity*length/viscosity = 10, safely laminar

        LidDrivenCavitySolver2D solver(n, n, length, length, viscosity, lidVelocity);
        const double dt = solver.stableTimeStep();
        const int steps = 800;
        for (int s = 0; s < steps; ++s) {
            solver.step(dt);
        }

        std::printf("Aether Workspace [cavity] - cavidade com tampa deslizante (Re=%.0f)\n",
                    lidVelocity * length / viscosity);
        std::printf("  grade %zux%zu, passos: %d, divergencia max: %.4f\n", n, n, steps,
                    solver.maxDivergence());

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

        const GLuint vertexShader = gl33::compileShader(gl33::kGlVertexShader, kVertexShaderSource);
        const GLuint fragmentShader = gl33::compileShader(gl33::kGlFragmentShader, kFragmentShaderSource);
        program_ = gl33::linkProgram(vertexShader, fragmentShader);
        gl33::glDeleteShader(vertexShader);
        gl33::glDeleteShader(fragmentShader);
        projectionLoc_ = gl33::glGetUniformLocation(program_, "uProjection");

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
        arrowFieldVertexCount_ = static_cast<GLsizei>(lineVertexData.size() / 5);

        for (const auto& path : streamlines) {
            streamlineFirst_.push_back(static_cast<GLint>(lineVertexData.size() / 5));
            for (const auto& p : path) {
                pushVertex(lineVertexData, p.x, p.y, 1.0f, 0.9f, 0.2f);
            }
            streamlineCount_.push_back(static_cast<GLsizei>(path.size()));
        }

        gl33::glGenVertexArrays(1, &vao_);
        gl33::glBindVertexArray(vao_);
        gl33::glGenBuffers(1, &vbo_);
        gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo_);
        gl33::glBufferData(gl33::kGlArrayBuffer,
                            static_cast<gl33::GLsizeiptr>(lineVertexData.size() * sizeof(float)),
                            lineVertexData.data(), gl33::kGlStaticDraw);
        const GLsizei stride = 5 * sizeof(float);
        gl33::glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
        gl33::glEnableVertexAttribArray(0);
        gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(2 * sizeof(float)));
        gl33::glEnableVertexAttribArray(1);
        gl33::glBindVertexArray(0);

        const double margin = 0.05;
        projection_ = Matrix4x4::ortho(-margin, length + margin, -margin, length + margin, -1.0, 1.0);
        return true;
    }

    void shutdown() override {
        gl33::glDeleteBuffers(1, &vbo_);
        gl33::glDeleteVertexArrays(1, &vao_);
        gl33::glDeleteProgram(program_);
        streamlineFirst_.clear();
        streamlineCount_.clear();
    }

    void handleInput(const aether::app::UiInput&, bool, double) override {}
    void update() override {}

    void renderScene(int, int) override {
        glClearColor(0.05f, 0.05f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        gl33::glUseProgram(program_);
        gl33::glUniformMatrix4fv(projectionLoc_, 1, GL_FALSE, projection_.data());
        gl33::glBindVertexArray(vao_);
        glDrawArrays(GL_LINES, 0, arrowFieldVertexCount_);
        for (std::size_t s = 0; s < streamlineFirst_.size(); ++s) {
            glDrawArrays(GL_LINE_STRIP, streamlineFirst_[s], streamlineCount_[s]);
        }
        gl33::glBindVertexArray(0);
    }

    void renderPanel(aether::app::Ui&) override {}

private:
    static constexpr const char* kVertexShaderSource = R"GLSL(
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

    static constexpr const char* kFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
}
)GLSL";

    GLuint program_ = 0;
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLint projectionLoc_ = -1;
    GLsizei arrowFieldVertexCount_ = 0;
    std::vector<GLint> streamlineFirst_;
    std::vector<GLsizei> streamlineCount_;
    aether::core::Matrix4x4 projection_;
};

} // namespace aether_app_modes

// ===========================================================================
// Mode: cavity3d -- the first viewer for any of this project's 3D solvers.
// StaggeredLidDrivenCavitySolver3D's velocity field as a 3D arrow field,
// colored by speed, on the mesh mode's perspective/orbit-camera skeleton.
// ===========================================================================

namespace aether_app_modes {

class Cavity3DMode : public Mode {
public:
    const char* label() const override { return "CAVIDADE 3D (cavity3d)"; }

    bool init() override {
        using aether::core::Matrix4x4;
        using aether::core::Vector3;
        using aether::solver::StaggeredLidDrivenCavitySolver3D;

        const std::size_t n = 10;
        const double length = 1.0;
        const double viscosity = 0.1;
        const double lidVelocity = 1.0;

        StaggeredLidDrivenCavitySolver3D solver(n, n, n, length, length, length, viscosity, lidVelocity);
        const double dt = solver.stableTimeStep();
        const int steps = 400;
        for (int s = 0; s < steps; ++s) {
            solver.step(dt);
        }

        std::printf("Aether Workspace [cavity3d] - cavidade 3D com tampa deslizante (Re=%.0f)\n",
                    lidVelocity * length / viscosity);
        std::printf("  grade %zux%zux%zu, passos: %d, divergencia max: %.6f\n", n, n, n, steps,
                    solver.maxDivergence());

        camera_.center = Vector3(length * 0.5, length * 0.5, length * 0.5);
        camera_.distance = length * 2.5;

        const GLuint vertexShader = gl33::compileShader(gl33::kGlVertexShader, kVertexShaderSource);
        const GLuint fragmentShader = gl33::compileShader(gl33::kGlFragmentShader, kFragmentShaderSource);
        program_ = gl33::linkProgram(vertexShader, fragmentShader);
        gl33::glDeleteShader(vertexShader);
        gl33::glDeleteShader(fragmentShader);
        mvpLoc_ = gl33::glGetUniformLocation(program_, "uMvp");

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
        boxVertexCount_ = static_cast<GLsizei>(vertexData.size() / 6);

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
        arrowVertexCount_ = static_cast<GLsizei>(vertexData.size() / 6 - boxVertexCount_);

        gl33::glGenVertexArrays(1, &vao_);
        gl33::glBindVertexArray(vao_);
        gl33::glGenBuffers(1, &vbo_);
        gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo_);
        gl33::glBufferData(gl33::kGlArrayBuffer,
                            static_cast<gl33::GLsizeiptr>(vertexData.size() * sizeof(float)), vertexData.data(),
                            gl33::kGlStaticDraw);
        const GLsizei stride = 6 * sizeof(float);
        gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
        gl33::glEnableVertexAttribArray(0);
        gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(3 * sizeof(float)));
        gl33::glEnableVertexAttribArray(1);
        gl33::glBindVertexArray(0);
        return true;
    }

    void shutdown() override {
        gl33::glDeleteBuffers(1, &vbo_);
        gl33::glDeleteVertexArrays(1, &vao_);
        gl33::glDeleteProgram(program_);
        dragging_ = false;
    }

    void handleInput(const aether::app::UiInput& input, bool uiWantsMouse, double wheelDelta) override {
        if (wheelDelta != 0.0) {
            const double factor = wheelDelta > 0 ? 0.9 : 1.1;
            camera_.distance = std::max(1e-3, camera_.distance * factor);
        }
        if (input.mousePressed && !uiWantsMouse) {
            dragging_ = true;
            lastMouseX_ = input.mouseX;
            lastMouseY_ = input.mouseY;
        }
        if (input.mouseReleased) {
            dragging_ = false;
        }
        if (dragging_ && input.mouseDown) {
            camera_.yaw += (input.mouseX - lastMouseX_) * 0.01;
            camera_.pitch = std::clamp(camera_.pitch + (input.mouseY - lastMouseY_) * 0.01, -1.5, 1.5);
            lastMouseX_ = input.mouseX;
            lastMouseY_ = input.mouseY;
        }
    }

    void update() override {}

    void renderScene(int width, int height) override {
        using aether::core::Matrix4x4;
        using aether::core::Vector3;

        glEnable(GL_DEPTH_TEST);
        glClearColor(0.08f, 0.08f, 0.1f, 1.0f);
        glLineWidth(1.5f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const double aspect = height > 0 ? static_cast<double>(width) / height : 1.0;
        const double nearPlane = std::max(0.01 * camera_.distance, 1e-3);
        const double farPlane = camera_.distance * 100.0 + 10.0;
        const Matrix4x4 projection = Matrix4x4::perspective(45.0 * 3.14159265358979323846 / 180.0, aspect,
                                                              nearPlane, farPlane);
        const Vector3 eye = camera_.eye();
        const Matrix4x4 view = Matrix4x4::lookAt(eye, camera_.center, Vector3(0.0, 1.0, 0.0));
        const Matrix4x4 mvp = projection * view;

        gl33::glUseProgram(program_);
        gl33::glUniformMatrix4fv(mvpLoc_, 1, GL_FALSE, mvp.data());

        gl33::glBindVertexArray(vao_);
        glDrawArrays(GL_LINES, 0, boxVertexCount_);
        glDrawArrays(GL_LINES, boxVertexCount_, arrowVertexCount_);
        gl33::glBindVertexArray(0);
    }

    void renderPanel(aether::app::Ui&) override {}

private:
    struct OrbitCamera {
        aether::core::Vector3 center;
        double distance = 3.0;
        double yaw = 0.7;
        double pitch = 0.5;

        aether::core::Vector3 eye() const {
            const double cp = std::cos(pitch);
            return center + aether::core::Vector3(distance * cp * std::cos(yaw), distance * std::sin(pitch),
                                                    distance * cp * std::sin(yaw));
        }
    };

    static constexpr const char* kVertexShaderSource = R"GLSL(
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

    static constexpr const char* kFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
}
)GLSL";

    OrbitCamera camera_;
    bool dragging_ = false;
    int lastMouseX_ = 0;
    int lastMouseY_ = 0;

    GLuint program_ = 0;
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLint mvpLoc_ = -1;
    GLsizei boxVertexCount_ = 0;
    GLsizei arrowVertexCount_ = 0;
};

} // namespace aether_app_modes

// ===========================================================================
// Mode: isosurface -- eddy-viscosity iso-surface from the 3D k-omega SST
// cavity, extracted via marchingCubes3D and rendered as a real shaded
// triangle mesh (mesh mode's lit pipeline), not arrows.
// ===========================================================================

namespace aether_app_modes {

class IsosurfaceMode : public Mode {
public:
    const char* label() const override { return "ISO-SUPERFICIE (isosurface)"; }

    bool init() override {
        using aether::core::Matrix4x4;
        using aether::core::Vector3;
        using aether::postprocessing::marchingCubes3D;
        using aether::postprocessing::Triangle3D;
        using aether::solver::KOmegaSSTLidDrivenCavitySolver3D;

        const std::size_t n = 12;
        const double length = 1.0;
        const double viscosity = 0.01;
        const double lidVelocity = 1.0;

        KOmegaSSTLidDrivenCavitySolver3D solver(n, n, n, length, length, length, viscosity, lidVelocity);
        double dt = solver.stableTimeStep();
        const int steps = 300;
        for (int s = 0; s < steps; ++s) {
            solver.step(dt);
            if (s % 100 == 0) {
                dt = solver.stableTimeStep();
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

        const double isoValue = 0.3 * std::max(maxNut, 1e-12);
        const auto triangles = marchingCubes3D(n, n, n, length, length, length, nut, isoValue);

        std::printf("Aether Workspace [isosurface] - iso-superficie de nu_t (k-omega SST 3D, Re=%.0f)\n",
                    lidVelocity * length / viscosity);
        std::printf("  grade %zux%zux%zu, passos: %d, max(nu_t)=%.6e, iso=%.6e, triangulos: %zu\n", n, n, n,
                    steps, maxNut, isoValue, triangles.size());

        camera_.center = Vector3(length * 0.5, length * 0.5, length * 0.5);
        camera_.distance = length * 2.5;

        const GLuint surfaceVs = gl33::compileShader(gl33::kGlVertexShader, kSurfaceVertexShaderSource);
        const GLuint surfaceFs = gl33::compileShader(gl33::kGlFragmentShader, kSurfaceFragmentShaderSource);
        surfaceProgram_ = gl33::linkProgram(surfaceVs, surfaceFs);
        gl33::glDeleteShader(surfaceVs);
        gl33::glDeleteShader(surfaceFs);
        surfaceMvpLoc_ = gl33::glGetUniformLocation(surfaceProgram_, "uMvp");
        lightDirLoc_ = gl33::glGetUniformLocation(surfaceProgram_, "uLightDir");
        baseColorLoc_ = gl33::glGetUniformLocation(surfaceProgram_, "uBaseColor");

        const GLuint lineVs = gl33::compileShader(gl33::kGlVertexShader, kLineVertexShaderSource);
        const GLuint lineFs = gl33::compileShader(gl33::kGlFragmentShader, kLineFragmentShaderSource);
        lineProgram_ = gl33::linkProgram(lineVs, lineFs);
        gl33::glDeleteShader(lineVs);
        gl33::glDeleteShader(lineFs);
        lineMvpLoc_ = gl33::glGetUniformLocation(lineProgram_, "uMvp");

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
        surfaceVertexCount_ = static_cast<GLsizei>(triangles.size() * 3);

        gl33::glGenVertexArrays(1, &surfaceVao_);
        gl33::glBindVertexArray(surfaceVao_);
        gl33::glGenBuffers(1, &surfaceVbo_);
        gl33::glBindBuffer(gl33::kGlArrayBuffer, surfaceVbo_);
        gl33::glBufferData(gl33::kGlArrayBuffer,
                            static_cast<gl33::GLsizeiptr>(surfaceVertexData.size() * sizeof(float)),
                            surfaceVertexData.empty() ? nullptr : surfaceVertexData.data(),
                            gl33::kGlStaticDraw);
        const GLsizei surfaceStride = 6 * sizeof(float);
        gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, surfaceStride, reinterpret_cast<void*>(0));
        gl33::glEnableVertexAttribArray(0);
        gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, surfaceStride,
                                     reinterpret_cast<void*>(3 * sizeof(float)));
        gl33::glEnableVertexAttribArray(1);
        gl33::glBindVertexArray(0);

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
        boxVertexCount_ = static_cast<GLsizei>(lineVertexData.size() / 6);

        gl33::glGenVertexArrays(1, &lineVao_);
        gl33::glBindVertexArray(lineVao_);
        gl33::glGenBuffers(1, &lineVbo_);
        gl33::glBindBuffer(gl33::kGlArrayBuffer, lineVbo_);
        gl33::glBufferData(gl33::kGlArrayBuffer,
                            static_cast<gl33::GLsizeiptr>(lineVertexData.size() * sizeof(float)),
                            lineVertexData.data(), gl33::kGlStaticDraw);
        const GLsizei lineStride = 6 * sizeof(float);
        gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, lineStride, reinterpret_cast<void*>(0));
        gl33::glEnableVertexAttribArray(0);
        gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, lineStride, reinterpret_cast<void*>(3 * sizeof(float)));
        gl33::glEnableVertexAttribArray(1);
        gl33::glBindVertexArray(0);
        return true;
    }

    void shutdown() override {
        gl33::glDeleteBuffers(1, &surfaceVbo_);
        gl33::glDeleteVertexArrays(1, &surfaceVao_);
        gl33::glDeleteProgram(surfaceProgram_);
        gl33::glDeleteBuffers(1, &lineVbo_);
        gl33::glDeleteVertexArrays(1, &lineVao_);
        gl33::glDeleteProgram(lineProgram_);
        dragging_ = false;
    }

    void handleInput(const aether::app::UiInput& input, bool uiWantsMouse, double wheelDelta) override {
        if (wheelDelta != 0.0) {
            const double factor = wheelDelta > 0 ? 0.9 : 1.1;
            camera_.distance = std::max(1e-3, camera_.distance * factor);
        }
        if (input.mousePressed && !uiWantsMouse) {
            dragging_ = true;
            lastMouseX_ = input.mouseX;
            lastMouseY_ = input.mouseY;
        }
        if (input.mouseReleased) {
            dragging_ = false;
        }
        if (dragging_ && input.mouseDown) {
            camera_.yaw += (input.mouseX - lastMouseX_) * 0.01;
            camera_.pitch = std::clamp(camera_.pitch + (input.mouseY - lastMouseY_) * 0.01, -1.5, 1.5);
            lastMouseX_ = input.mouseX;
            lastMouseY_ = input.mouseY;
        }
    }

    void update() override {}

    void renderScene(int width, int height) override {
        using aether::core::Matrix4x4;
        using aether::core::Vector3;

        glEnable(GL_DEPTH_TEST);
        glClearColor(0.08f, 0.08f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const double aspect = height > 0 ? static_cast<double>(width) / height : 1.0;
        const double nearPlane = std::max(0.01 * camera_.distance, 1e-3);
        const double farPlane = camera_.distance * 100.0 + 10.0;
        const Matrix4x4 projection = Matrix4x4::perspective(45.0 * 3.14159265358979323846 / 180.0, aspect,
                                                              nearPlane, farPlane);
        const Vector3 eye = camera_.eye();
        const Matrix4x4 view = Matrix4x4::lookAt(eye, camera_.center, Vector3(0.0, 1.0, 0.0));
        const Matrix4x4 mvp = projection * view;

        if (surfaceVertexCount_ > 0) {
            gl33::glUseProgram(surfaceProgram_);
            gl33::glUniformMatrix4fv(surfaceMvpLoc_, 1, GL_FALSE, mvp.data());
            gl33::glUniform3f(lightDirLoc_, 1.0f, 1.0f, 1.0f);
            gl33::glUniform3f(baseColorLoc_, 0.85f, 0.55f, 0.2f);
            gl33::glBindVertexArray(surfaceVao_);
            glDrawArrays(GL_TRIANGLES, 0, surfaceVertexCount_);
            gl33::glBindVertexArray(0);
        }

        gl33::glUseProgram(lineProgram_);
        gl33::glUniformMatrix4fv(lineMvpLoc_, 1, GL_FALSE, mvp.data());
        gl33::glBindVertexArray(lineVao_);
        glDrawArrays(GL_LINES, 0, boxVertexCount_);
        gl33::glBindVertexArray(0);
    }

    void renderPanel(aether::app::Ui&) override {}

private:
    struct OrbitCamera {
        aether::core::Vector3 center;
        double distance = 3.0;
        double yaw = 0.7;
        double pitch = 0.5;

        aether::core::Vector3 eye() const {
            const double cp = std::cos(pitch);
            return center + aether::core::Vector3(distance * cp * std::cos(yaw), distance * std::sin(pitch),
                                                    distance * cp * std::sin(yaw));
        }
    };

    static constexpr const char* kSurfaceVertexShaderSource = R"GLSL(
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

    static constexpr const char* kSurfaceFragmentShaderSource = R"GLSL(
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

    static constexpr const char* kLineVertexShaderSource = R"GLSL(
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

    static constexpr const char* kLineFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
}
)GLSL";

    OrbitCamera camera_;
    bool dragging_ = false;
    int lastMouseX_ = 0;
    int lastMouseY_ = 0;

    GLuint surfaceProgram_ = 0;
    GLuint surfaceVao_ = 0;
    GLuint surfaceVbo_ = 0;
    GLint surfaceMvpLoc_ = -1;
    GLint lightDirLoc_ = -1;
    GLint baseColorLoc_ = -1;
    GLsizei surfaceVertexCount_ = 0;

    GLuint lineProgram_ = 0;
    GLuint lineVao_ = 0;
    GLuint lineVbo_ = 0;
    GLint lineMvpLoc_ = -1;
    GLsizei boxVertexCount_ = 0;
};

} // namespace aether_app_modes

// ===========================================================================
// Mode: turbulence (formerly apps/turbulence_viewer) -- u+ vs ln(y+).
// ===========================================================================

namespace aether_app_modes {

class TurbulenceMode : public Mode {
public:
    const char* label() const override { return "TURBULENCIA (turbulence)"; }

    bool init() override {
        using aether::core::Matrix4x4;
        using aether::solver::KEpsilonChannelFlowSolver1D;
        using aether::solver::KOmegaSSTChannelFlowSolver1D;
        using aether::solver::MixingLengthChannelFlowSolver1D;

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
        std::printf("Aether Workspace [turbulence] - u+ vs ln(y+), Re_tau ~ %.0f\n", reTau);
        std::printf("  comprimento de mistura: %zu iteracoes, u_tau=%.5f\n", mixingLengthIterations,
                    uTauMixingLength);
        std::printf("  k-epsilon:              %zu iteracoes, u_tau=%.5f\n", kEpsilonIterations, uTauKEpsilon);
        std::printf("  k-omega SST:            %zu iteracoes, u_tau=%.5f\n", kOmegaSSTIterations,
                    uTauKOmegaSST);

        const double xMin = logYPlus(1.0);
        const double xMax = logYPlus(2000.0);
        const double yMin = 0.0;
        const double yMax = 26.0;

        const GLuint vertexShader = gl33::compileShader(gl33::kGlVertexShader, kVertexShaderSource);
        const GLuint fragmentShader = gl33::compileShader(gl33::kGlFragmentShader, kFragmentShaderSource);
        program_ = gl33::linkProgram(vertexShader, fragmentShader);
        gl33::glDeleteShader(vertexShader);
        gl33::glDeleteShader(fragmentShader);
        projectionLoc_ = gl33::glGetUniformLocation(program_, "uProjection");

        std::vector<float> vertexData;
        ranges_.clear();
        auto beginRange = [&]() { return static_cast<GLint>(vertexData.size() / 5); };
        auto endRange = [&](GLenum primitive, GLint first) {
            ranges_.push_back({primitive, first, static_cast<GLsizei>(vertexData.size() / 5 - first)});
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

        gl33::glGenVertexArrays(1, &vao_);
        gl33::glBindVertexArray(vao_);
        gl33::glGenBuffers(1, &vbo_);
        gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo_);
        gl33::glBufferData(gl33::kGlArrayBuffer,
                            static_cast<gl33::GLsizeiptr>(vertexData.size() * sizeof(float)), vertexData.data(),
                            gl33::kGlStaticDraw);
        const GLsizei stride = 5 * sizeof(float);
        gl33::glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
        gl33::glEnableVertexAttribArray(0);
        gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(2 * sizeof(float)));
        gl33::glEnableVertexAttribArray(1);
        gl33::glBindVertexArray(0);

        const double marginFrac = 0.08;
        const double xSpan = xMax - xMin;
        const double ySpan = yMax - yMin;
        projection_ = Matrix4x4::ortho(xMin - marginFrac * xSpan, xMax + marginFrac * xSpan,
                                        yMin - marginFrac * ySpan, yMax + marginFrac * ySpan, -1.0, 1.0);
        return true;
    }

    void shutdown() override {
        gl33::glDeleteBuffers(1, &vbo_);
        gl33::glDeleteVertexArrays(1, &vao_);
        gl33::glDeleteProgram(program_);
        ranges_.clear();
    }

    void handleInput(const aether::app::UiInput&, bool, double) override {}
    void update() override {}

    void renderScene(int, int) override {
        glClearColor(0.08f, 0.08f, 0.1f, 1.0f);
        glLineWidth(2.0f);
        glPointSize(6.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        gl33::glUseProgram(program_);
        gl33::glUniformMatrix4fv(projectionLoc_, 1, GL_FALSE, projection_.data());
        gl33::glBindVertexArray(vao_);
        for (const DrawRange& range : ranges_) {
            glDrawArrays(range.primitive, range.first, range.count);
        }
        gl33::glBindVertexArray(0);
    }

    void renderPanel(aether::app::Ui&) override {}

private:
    struct DrawRange {
        GLenum primitive;
        GLint first;
        GLsizei count;
    };

    static double logYPlus(double yPlus) { return std::log(std::max(yPlus, 1e-6)); }

    static constexpr const char* kVertexShaderSource = R"GLSL(
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

    static constexpr const char* kFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
}
)GLSL";

    static constexpr double kKarman = 0.41;
    static constexpr double kLogLawB = 5.0;

    GLuint program_ = 0;
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLint projectionLoc_ = -1;
    std::vector<DrawRange> ranges_;
    aether::core::Matrix4x4 projection_;
};

} // namespace aether_app_modes

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

namespace aether_app_modes {

class SimMode : public Mode {
public:
    const char* label() const override { return "SIMULACAO 2D (sim)"; }

    bool init() override {
        std::printf("Aether Workspace - modo 'sim' (Modulo 9: UI com paineis)\n");
        std::printf("  painel a esquerda: fechamento, parametros, run/pause/passo/reiniciar\n");
        std::printf("  o campo e' redesenhado a cada quadro a partir do solver ao vivo\n\n");
        std::fflush(stdout);

        sim_ = Simulation{};
        sim_.rebuild();

        running_ = true;
        showPressure_ = false;
        resolution_ = static_cast<double>(sim_.n);
        viscosity_ = sim_.viscosity;
        lidVelocity_ = sim_.lidVelocity;
        return true;
    }

    void shutdown() override {
        sim_.laminar.reset();
        sim_.mixingLength.reset();
        sim_.kEpsilon.reset();
        sim_.kOmegaSST.reset();
    }

    void handleInput(const aether::app::UiInput&, bool, double) override {}

    void update() override {
        if (running_) {
            for (int s = 0; s < kStepsPerFrame; ++s) {
                sim_.step();
            }
        }
    }

    void renderScene(int width, int height) override {
        width_ = width;
        height_ = height;
        glClearColor(0.06f, 0.07f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }

    void renderPanel(aether::app::Ui& ui) override {
        const int width = width_;
        const int height = height_;

        const int panelWidth = 260;
        const int margin = kSidebarMargin;
        const int panelRight = kModePanelX + panelWidth;
        const int available = std::min(width - panelRight - 2 * margin, height - 2 * margin);
        const int fieldSize = std::max(available, 64);
        const int fieldX = panelRight + margin;
        const int fieldY = (height - fieldSize) / 2;

        double maxValue = 1e-12;
        double minValue = 0.0;
        for (std::size_t j = 0; j < sim_.n; ++j) {
            for (std::size_t i = 0; i < sim_.n; ++i) {
                const double value = showPressure_ ? sim_.pressureAt(i, j) : sim_.speedAt(i, j);
                maxValue = std::max(maxValue, value);
                minValue = std::min(minValue, value);
            }
        }
        const double span = (maxValue - minValue) > 1e-12 ? (maxValue - minValue) : 1.0;

        const int cell = std::max(fieldSize / static_cast<int>(sim_.n), 1);
        for (std::size_t j = 0; j < sim_.n; ++j) {
            for (std::size_t i = 0; i < sim_.n; ++i) {
                const double value = showPressure_ ? sim_.pressureAt(i, j) : sim_.speedAt(i, j);
                float r = 0.0f;
                float g = 0.0f;
                float b = 0.0f;
                colorForValue((value - minValue) / span, r, g, b);
                const int px = fieldX + static_cast<int>(i) * cell;
                const int py = fieldY + (static_cast<int>(sim_.n) - 1 - static_cast<int>(j)) * cell;
                ui.drawRect(px, py, cell, cell, aether::app::UiColor{r, g, b, 1.0f});
            }
        }
        ui.drawText(fieldX, fieldY - 16,
                    showPressure_ ? "pressao (azul=baixa, vermelho=alta)"
                                  : "velocidade |u| (azul=baixa, vermelho=alta)",
                    aether::app::UiColor{0.60f, 0.63f, 0.67f, 1.0f});

        ui.beginPanel(kModePanelX, margin, panelWidth, "AETHER - CAVIDADE 2D");

        ui.label("Fechamento de turbulencia");
        for (int c = 0; c < 4; ++c) {
            const auto candidate = static_cast<Closure>(c);
            std::string caption = (sim_.closure == candidate) ? "[x] " : "[ ] ";
            caption += closureName(candidate);
            if (ui.button(caption) && sim_.closure != candidate) {
                sim_.closure = candidate;
                sim_.rebuild();
                running_ = false;
            }
        }

        ui.separator();
        ui.label("Parametros");
        bool needsRebuild = false;
        needsRebuild |= ui.slider("resolucao", &resolution_, 16, 80);
        needsRebuild |= ui.textField("viscosidade", &viscosity_, 0.0005, 0.5);
        needsRebuild |= ui.textField("veloc. da tampa", &lidVelocity_, 0.0, 5.0);
        if (needsRebuild) {
            sim_.n = static_cast<std::size_t>(resolution_);
            sim_.viscosity = viscosity_;
            sim_.lidVelocity = lidVelocity_;
            sim_.rebuild();
            running_ = false;
        }

        ui.separator();
        ui.label("Simulacao");
        if (ui.button(running_ ? "PAUSAR" : "RODAR")) {
            running_ = !running_;
        }
        if (ui.button("+1 PASSO")) {
            sim_.step();
        }
        if (ui.button("REINICIAR")) {
            sim_.rebuild();
            running_ = false;
        }
        ui.checkbox("mostrar pressao", &showPressure_);

        ui.separator();
        ui.label("Diagnostico");
        char buffer[128];
        std::snprintf(buffer, sizeof(buffer), "Re = %.0f", sim_.lidVelocity / sim_.viscosity);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "malha: %zux%zu", sim_.n, sim_.n);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "passos: %lld", sim_.steps);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "t = %.4f s", sim_.time());
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "div max = %.2e", sim_.maxDivergence());
        ui.label(buffer);

        ui.endPanel();
    }

private:
    enum class Closure { Laminar = 0, MixingLength, KEpsilon, KOmegaSST };

    static const char* closureName(Closure closure) {
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
            const auto advance = [](auto& solver) { solver.step(solver.stableTimeStep()); };
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

    static constexpr int kStepsPerFrame = 4;

    Simulation sim_;
    bool running_ = false;
    bool showPressure_ = false;
    double resolution_ = 40.0;
    double viscosity_ = 0.01;
    double lidVelocity_ = 1.0;
    int width_ = 1280;
    int height_ = 800;
};

} // namespace aether_app_modes


// ===========================================================================
// Mode: sim3d (Module 9.5) -- the sim mode's control panel over the real 3D
// staggered cavity, across all six closures this project has (laminar,
// mixing-length, k-epsilon, k-omega SST, Smagorinsky LES, SST-DES).
//
// **Why this needs a worker thread and sim_mode's 2D panel did not**: a 3D
// step costs O(n^3) instead of O(n^2), and the k-omega/DES/k-epsilon
// constructors each run a 400-step mixing-length primer internally before
// the panel can show anything -- both cheap enough at 2D resolutions to sit
// directly in the render loop, neither cheap enough in 3D to keep a 60fps
// UI responsive. So the solver now lives on its own thread: the worker
// advances it (or rebuilds it) continuously, guarded by a mutex, and the
// render thread only ever takes a brief try_lock() once per frame to copy
// out whatever it needs to draw. If the worker is mid-step or mid-rebuild
// when a frame wants to render, try_lock() fails, and the frame simply
// redraws the previous snapshot -- the panel and the orbit camera stay
// responsive even while a rebuild (which can take a second or two at these
// resolutions) is in flight on the other thread. This is the first
// multi-threaded code in the project.
// ===========================================================================

namespace aether_app_modes {

class Sim3DMode : public Mode {
public:
    const char* label() const override { return "SIMULACAO 3D (sim3d)"; }

    bool init() override {
        using aether::core::Vector3;

        std::printf("Aether Workspace - modo 'sim3d' (Modulo 9.5: painel sobre a cavidade 3D)\n");
        std::printf("  o solver 3D roda numa thread propria; arraste fora do painel gira a camera\n\n");
        std::fflush(stdout);

        const GLuint vertexShader = gl33::compileShader(gl33::kGlVertexShader, kVertexShaderSource);
        const GLuint fragmentShader = gl33::compileShader(gl33::kGlFragmentShader, kFragmentShaderSource);
        program_ = gl33::linkProgram(vertexShader, fragmentShader);
        gl33::glDeleteShader(vertexShader);
        gl33::glDeleteShader(fragmentShader);
        mvpLoc_ = gl33::glGetUniformLocation(program_, "uMvp");

        gl33::glGenVertexArrays(1, &vao_);
        gl33::glBindVertexArray(vao_);
        gl33::glGenBuffers(1, &vbo_);
        gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo_);
        constexpr GLsizei kStride = 6 * sizeof(float);
        gl33::glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, kStride, reinterpret_cast<void*>(0));
        gl33::glEnableVertexAttribArray(0);
        gl33::glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, kStride,
                                     reinterpret_cast<void*>(3 * sizeof(float)));
        gl33::glEnableVertexAttribArray(1);
        gl33::glBindVertexArray(0);

        sim_ = std::make_unique<Simulation3D>();
        workerQuit_.store(false, std::memory_order_relaxed);
        worker_ = std::thread(workerMain, sim_.get(), &workerQuit_);

        closure_ = Closure3D::Laminar;
        resolution_ = 10.0;
        viscosity_ = 0.05;
        lidVelocity_ = 1.0;
        running_ = true;
        sim_->setRunning(running_);
        sim_->requestRebuild(closure_, static_cast<std::size_t>(resolution_), viscosity_, lidVelocity_);

        camera_ = OrbitCamera{};
        camera_.center = Vector3(0.5, 0.5, 0.5);
        camera_.distance = 2.6;
        dragging_ = false;

        snapshot_ = Snapshot3D{};
        return true;
    }

    void shutdown() override {
        workerQuit_.store(true, std::memory_order_relaxed);
        if (worker_.joinable()) {
            worker_.join();
        }
        sim_.reset();

        gl33::glDeleteBuffers(1, &vbo_);
        gl33::glDeleteVertexArrays(1, &vao_);
        gl33::glDeleteProgram(program_);
        dragging_ = false;
    }

    void handleInput(const aether::app::UiInput& input, bool uiWantsMouse, double wheelDelta) override {
        if (wheelDelta != 0.0) {
            const double factor = wheelDelta > 0 ? 0.9 : 1.1;
            camera_.distance = std::max(1e-3, camera_.distance * factor);
        }
        if (input.mousePressed && !uiWantsMouse) {
            dragging_ = true;
            lastMouseX_ = input.mouseX;
            lastMouseY_ = input.mouseY;
        }
        if (input.mouseReleased) {
            dragging_ = false;
        }
        if (dragging_ && input.mouseDown) {
            camera_.yaw += (input.mouseX - lastMouseX_) * 0.01;
            camera_.pitch = std::clamp(camera_.pitch + (input.mouseY - lastMouseY_) * 0.01, -1.5, 1.5);
            lastMouseX_ = input.mouseX;
            lastMouseY_ = input.mouseY;
        }
    }

    void update() override { sim_->trySnapshot(snapshot_); }

    void renderScene(int width, int height) override {
        using aether::core::Matrix4x4;
        using aether::core::Vector3;

        glEnable(GL_DEPTH_TEST);
        glLineWidth(1.5f);
        glClearColor(0.06f, 0.07f, 0.08f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        std::vector<float> vertexData;
        const auto pushVertex3 = [&vertexData](double x, double y, double z, float r, float g, float b) {
            vertexData.push_back(static_cast<float>(x));
            vertexData.push_back(static_cast<float>(y));
            vertexData.push_back(static_cast<float>(z));
            vertexData.push_back(r);
            vertexData.push_back(g);
            vertexData.push_back(b);
        };
        {
            const float c[3] = {0.5f, 0.5f, 0.55f};
            const double corners[8][3] = {
                {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
            };
            const int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                                       {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
            for (const auto& e : edges) {
                pushVertex3(corners[e[0]][0], corners[e[0]][1], corners[e[0]][2], c[0], c[1], c[2]);
                pushVertex3(corners[e[1]][0], corners[e[1]][1], corners[e[1]][2], c[0], c[1], c[2]);
            }
        }
        const auto boxVertexCount = static_cast<GLsizei>(vertexData.size() / 6);

        if (snapshot_.n > 0) {
            const std::size_t n = snapshot_.n;
            const double cellSize = 1.0 / static_cast<double>(n);
            double maxSpeed = 1e-9;
            for (std::size_t idx = 0; idx < n * n * n; ++idx) {
                const double speed = std::sqrt(snapshot_.u[idx] * snapshot_.u[idx] +
                                                snapshot_.v[idx] * snapshot_.v[idx] +
                                                snapshot_.w[idx] * snapshot_.w[idx]);
                maxSpeed = std::max(maxSpeed, speed);
            }
            const double arrowScale = cellSize * 0.9 / maxSpeed;
            std::size_t idx = 0;
            for (std::size_t k = 0; k < n; ++k) {
                for (std::size_t j = 0; j < n; ++j) {
                    for (std::size_t i = 0; i < n; ++i, ++idx) {
                        const double cx = (static_cast<double>(i) + 0.5) * cellSize;
                        const double cy = (static_cast<double>(j) + 0.5) * cellSize;
                        const double cz = (static_cast<double>(k) + 0.5) * cellSize;
                        const double uc = snapshot_.u[idx];
                        const double vc = snapshot_.v[idx];
                        const double wc = snapshot_.w[idx];
                        const double speed = std::sqrt(uc * uc + vc * vc + wc * wc);
                        float r = 0.0f;
                        float g = 0.0f;
                        float b = 0.0f;
                        colorForValue(speed / maxSpeed, r, g, b);
                        pushVertex3(cx, cy, cz, r, g, b);
                        pushVertex3(cx + uc * arrowScale, cy + vc * arrowScale, cz + wc * arrowScale, r, g, b);
                    }
                }
            }
        }
        const auto arrowVertexCount = static_cast<GLsizei>(vertexData.size() / 6 - boxVertexCount);

        gl33::glBindBuffer(gl33::kGlArrayBuffer, vbo_);
        gl33::glBufferData(gl33::kGlArrayBuffer,
                            static_cast<gl33::GLsizeiptr>(vertexData.size() * sizeof(float)), vertexData.data(),
                            gl33::kGlStaticDraw);

        const double aspect = height > 0 ? static_cast<double>(width) / height : 1.0;
        const double nearPlane = std::max(0.01 * camera_.distance, 1e-3);
        const double farPlane = camera_.distance * 100.0 + 10.0;
        const Matrix4x4 projection = Matrix4x4::perspective(45.0 * 3.14159265358979323846 / 180.0, aspect,
                                                              nearPlane, farPlane);
        const Matrix4x4 view = Matrix4x4::lookAt(camera_.eye(), camera_.center, Vector3(0.0, 1.0, 0.0));
        const Matrix4x4 mvp = projection * view;

        gl33::glUseProgram(program_);
        gl33::glUniformMatrix4fv(mvpLoc_, 1, GL_FALSE, mvp.data());
        gl33::glBindVertexArray(vao_);
        glDrawArrays(GL_LINES, 0, boxVertexCount);
        glDrawArrays(GL_LINES, boxVertexCount, arrowVertexCount);
        gl33::glBindVertexArray(0);
    }

    void renderPanel(aether::app::Ui& ui) override {
        const int panelWidth = 260;
        const int margin = kSidebarMargin;
        ui.beginPanel(kModePanelX, margin, panelWidth, "AETHER - CAVIDADE 3D");

        ui.label("Fechamento de turbulencia");
        bool closureChanged = false;
        for (int c = 0; c <= static_cast<int>(Closure3D::DES); ++c) {
            const auto candidate = static_cast<Closure3D>(c);
            std::string caption = (closure_ == candidate) ? "[x] " : "[ ] ";
            caption += closureName3D(candidate);
            if (ui.button(caption) && closure_ != candidate) {
                closure_ = candidate;
                closureChanged = true;
            }
        }

        ui.separator();
        ui.label("Parametros");
        bool paramsChanged = false;
        paramsChanged |= ui.slider("resolucao", &resolution_, 6, 16);
        paramsChanged |= ui.textField("viscosidade", &viscosity_, 0.002, 0.5);
        paramsChanged |= ui.textField("veloc. da tampa", &lidVelocity_, 0.0, 5.0);
        if (closureChanged || paramsChanged) {
            sim_->requestRebuild(closure_, static_cast<std::size_t>(resolution_), viscosity_, lidVelocity_);
            running_ = false;
            sim_->setRunning(false);
        }

        ui.separator();
        ui.label("Simulacao");
        if (ui.button(running_ ? "PAUSAR" : "RODAR")) {
            running_ = !running_;
            sim_->setRunning(running_);
        }
        if (ui.button("+1 PASSO")) {
            sim_->requestStepOnce();
        }
        if (ui.button("REINICIAR")) {
            sim_->requestRebuild(closure_, static_cast<std::size_t>(resolution_), viscosity_, lidVelocity_);
            running_ = false;
            sim_->setRunning(false);
        }

        ui.separator();
        ui.label("Diagnostico");
        char buffer[128];
        std::snprintf(buffer, sizeof(buffer), "Re = %.0f", lidVelocity_ / std::max(viscosity_, 1e-9));
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "malha: %zux%zux%zu", snapshot_.n, snapshot_.n, snapshot_.n);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "passos: %lld", snapshot_.steps);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "t = %.4f s", snapshot_.time);
        ui.label(buffer);
        std::snprintf(buffer, sizeof(buffer), "div max = %.2e", snapshot_.maxDivergence);
        ui.label(buffer);
        ui.label(snapshot_.n > 0 ? "" : "(construindo...)");

        ui.endPanel();
    }

private:
    struct OrbitCamera {
        aether::core::Vector3 center;
        double distance = 3.0;
        double yaw = 0.7;
        double pitch = 0.5;

        aether::core::Vector3 eye() const {
            const double cp = std::cos(pitch);
            return center + aether::core::Vector3(distance * cp * std::cos(yaw), distance * std::sin(pitch),
                                                    distance * cp * std::sin(yaw));
        }
    };

    enum class Closure3D { Laminar = 0, MixingLength, KEpsilon, KOmegaSST, LES, DES };

    static const char* closureName3D(Closure3D closure) {
        switch (closure) {
        case Closure3D::Laminar:
            return "laminar";
        case Closure3D::MixingLength:
            return "mixing-length";
        case Closure3D::KEpsilon:
            return "k-epsilon";
        case Closure3D::KOmegaSST:
            return "k-omega SST";
        case Closure3D::LES:
            return "LES Smagorinsky";
        case Closure3D::DES:
            return "SST-DES";
        }
        return "?";
    }

    // One cell-centered velocity field, copied out under the lock so the
    // render thread never touches the live solver directly.
    struct Snapshot3D {
        std::size_t n = 0;
        std::vector<double> u;
        std::vector<double> v;
        std::vector<double> w;
        long long steps = 0;
        double time = 0.0;
        double maxDivergence = 0.0;
    };

    // Owns exactly one of the six 3D closures at a time, wrapped for
    // cross-thread access. Every method that touches the solver pointers or
    // the field arrays takes mutex_; the only exception is the plain
    // atomics (running_, pending flags), which don't need it.
    //
    // Threading contract: rebuildIfRequested() and stepOnce() are called
    // ONLY from the worker thread. requestRebuild(), requestStepOnce() and
    // setRunning() are called ONLY from the render/UI thread. trySnapshot()
    // is called ONLY from the render thread and never blocks -- it is the
    // one place a plain lock_guard would be wrong, since a blocking read
    // here would freeze the panel for as long as a rebuild takes.
    class Simulation3D {
    public:
        void requestRebuild(Closure3D closure, std::size_t n, double viscosity, double lidVelocity) {
            std::lock_guard<std::mutex> lock(mutex_);
            closure_ = closure;
            n_ = n;
            viscosity_ = viscosity;
            lidVelocity_ = lidVelocity;
            rebuildRequested_ = true;
        }
        void setRunning(bool running) { running_.store(running, std::memory_order_relaxed); }
        bool running() const { return running_.load(std::memory_order_relaxed); }
        void requestStepOnce() { pendingSingleStep_.store(true, std::memory_order_relaxed); }

        // --- worker-thread-only ---
        void rebuildIfRequested() {
            if (!rebuildRequested_.exchange(false)) {
                return;
            }
            std::lock_guard<std::mutex> lock(mutex_);
            ready_ = false;
            laminar_.reset();
            mixingLength_.reset();
            kEpsilon_.reset();
            kOmegaSST_.reset();
            les_.reset();
            des_.reset();
            constexpr double length = 1.0;
            switch (closure_) {
            case Closure3D::Laminar:
                laminar_ = std::make_unique<aether::solver::StaggeredLidDrivenCavitySolver3D>(
                    n_, n_, n_, length, length, length, viscosity_, lidVelocity_);
                break;
            case Closure3D::MixingLength:
                mixingLength_ = std::make_unique<aether::solver::MixingLengthLidDrivenCavitySolver3D>(
                    n_, n_, n_, length, length, length, viscosity_, lidVelocity_);
                break;
            case Closure3D::KEpsilon:
                kEpsilon_ = std::make_unique<aether::solver::KEpsilonLidDrivenCavitySolver3D>(
                    n_, n_, n_, length, length, length, viscosity_, lidVelocity_);
                break;
            case Closure3D::KOmegaSST:
                kOmegaSST_ = std::make_unique<aether::solver::KOmegaSSTLidDrivenCavitySolver3D>(
                    n_, n_, n_, length, length, length, viscosity_, lidVelocity_);
                break;
            case Closure3D::LES:
                les_ = std::make_unique<aether::solver::SmagorinskyLesLidDrivenCavitySolver3D>(
                    n_, n_, n_, length, length, length, viscosity_, lidVelocity_);
                break;
            case Closure3D::DES:
                des_ = std::make_unique<aether::solver::DesSstLidDrivenCavitySolver3D>(
                    n_, n_, n_, length, length, length, viscosity_, lidVelocity_);
                break;
            }
            steps_ = 0;
            time_ = 0.0;
            ready_ = true;
        }

        void stepOnce() {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!ready_) {
                return;
            }
            visitMutable([](auto& solver) { solver.step(solver.stableTimeStep()); });
            ++steps_;
        }

        bool consumeSingleStepRequest() {
            return pendingSingleStep_.exchange(false, std::memory_order_relaxed);
        }

        // --- render-thread-only ---
        // Never blocks: returns false (snapshot left untouched) if the
        // worker currently holds the lock, so a slow rebuild never stalls a
        // frame.
        bool trySnapshot(Snapshot3D& out) const {
            std::unique_lock<std::mutex> lock(mutex_, std::try_to_lock);
            if (!lock.owns_lock() || !ready_) {
                return false;
            }
            out.n = n_;
            const std::size_t total = n_ * n_ * n_;
            out.u.resize(total);
            out.v.resize(total);
            out.w.resize(total);
            std::size_t idx = 0;
            visit([&](const auto& solver) {
                for (std::size_t k = 0; k < n_; ++k) {
                    for (std::size_t j = 0; j < n_; ++j) {
                        for (std::size_t i = 0; i < n_; ++i, ++idx) {
                            out.u[idx] = 0.5 * (solver.u(i, j, k) + solver.u(i + 1, j, k));
                            out.v[idx] = 0.5 * (solver.v(i, j, k) + solver.v(i, j + 1, k));
                            out.w[idx] = 0.5 * (solver.w(i, j, k) + solver.w(i, j, k + 1));
                        }
                    }
                }
            });
            out.steps = steps_;
            out.time = visit([](const auto& solver) { return solver.time(); });
            out.maxDivergence = visit([](const auto& solver) { return solver.maxDivergence(); });
            return true;
        }

    private:
        template <typename Fn>
        auto visit(Fn&& fn) const
            -> decltype(fn(std::declval<aether::solver::StaggeredLidDrivenCavitySolver3D&>())) {
            if (laminar_) return fn(*laminar_);
            if (mixingLength_) return fn(*mixingLength_);
            if (kEpsilon_) return fn(*kEpsilon_);
            if (kOmegaSST_) return fn(*kOmegaSST_);
            if (les_) return fn(*les_);
            return fn(*des_);
        }
        template <typename Fn> void visitMutable(Fn&& fn) {
            if (laminar_) { fn(*laminar_); return; }
            if (mixingLength_) { fn(*mixingLength_); return; }
            if (kEpsilon_) { fn(*kEpsilon_); return; }
            if (kOmegaSST_) { fn(*kOmegaSST_); return; }
            if (les_) { fn(*les_); return; }
            fn(*des_);
        }

        mutable std::mutex mutex_;
        std::atomic<bool> running_{false};
        std::atomic<bool> rebuildRequested_{true};
        std::atomic<bool> pendingSingleStep_{false};
        bool ready_ = false;

        Closure3D closure_ = Closure3D::Laminar;
        std::size_t n_ = 10;
        double viscosity_ = 0.05;
        double lidVelocity_ = 1.0;
        long long steps_ = 0;
        double time_ = 0.0;

        std::unique_ptr<aether::solver::StaggeredLidDrivenCavitySolver3D> laminar_;
        std::unique_ptr<aether::solver::MixingLengthLidDrivenCavitySolver3D> mixingLength_;
        std::unique_ptr<aether::solver::KEpsilonLidDrivenCavitySolver3D> kEpsilon_;
        std::unique_ptr<aether::solver::KOmegaSSTLidDrivenCavitySolver3D> kOmegaSST_;
        std::unique_ptr<aether::solver::SmagorinskyLesLidDrivenCavitySolver3D> les_;
        std::unique_ptr<aether::solver::DesSstLidDrivenCavitySolver3D> des_;
    };

    static void workerMain(Simulation3D* sim, std::atomic<bool>* quit) {
        while (!quit->load(std::memory_order_relaxed)) {
            sim->rebuildIfRequested();
            if (sim->consumeSingleStepRequest() || sim->running()) {
                sim->stepOnce();
            } else {
                std::this_thread::sleep_for(std::chrono::milliseconds(5));
            }
        }
    }

    static constexpr const char* kVertexShaderSource = R"GLSL(
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

    static constexpr const char* kFragmentShaderSource = R"GLSL(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
}
)GLSL";

    std::unique_ptr<Simulation3D> sim_;
    std::thread worker_;
    std::atomic<bool> workerQuit_{false};

    Closure3D closure_ = Closure3D::Laminar;
    double resolution_ = 10.0;
    double viscosity_ = 0.05;
    double lidVelocity_ = 1.0;
    bool running_ = true;

    OrbitCamera camera_;
    bool dragging_ = false;
    int lastMouseX_ = 0;
    int lastMouseY_ = 0;

    Snapshot3D snapshot_;

    GLuint program_ = 0;
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLint mvpLoc_ = -1;
};

} // namespace aether_app_modes

namespace aether_app_modes {

namespace {
HWND g_workspaceWindow = nullptr;
int g_workspaceWidth = 1280;
int g_workspaceHeight = 800;
aether::app::UiInput g_workspaceInput;
bool g_workspacePendingPress = false;
bool g_workspacePendingRelease = false;
std::string g_workspacePendingText;
double g_workspacePendingWheelDelta = 0.0;

LRESULT CALLBACK workspaceWndProc(HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam) {
    switch (message) {
    case WM_DESTROY:
        PostQuitMessage(0);
        return 0;
    case WM_SIZE:
        g_workspaceWidth = LOWORD(lParam);
        g_workspaceHeight = HIWORD(lParam);
        if (wglGetCurrentContext()) {
            glViewport(0, 0, g_workspaceWidth, g_workspaceHeight);
        }
        return 0;
    case WM_MOUSEMOVE:
        g_workspaceInput.mouseX = GET_X_LPARAM(lParam);
        g_workspaceInput.mouseY = GET_Y_LPARAM(lParam);
        return 0;
    case WM_LBUTTONDOWN:
        g_workspaceInput.mouseDown = true;
        g_workspacePendingPress = true;
        SetCapture(hwnd);
        return 0;
    case WM_LBUTTONUP:
        g_workspaceInput.mouseDown = false;
        g_workspacePendingRelease = true;
        ReleaseCapture();
        return 0;
    case WM_MOUSEWHEEL: {
        const short delta = GET_WHEEL_DELTA_WPARAM(wParam);
        g_workspacePendingWheelDelta = delta > 0 ? 1.0 : -1.0;
        return 0;
    }
    case WM_CHAR:
        // TranslateMessage() (already called in this app's message loop)
        // turns WM_KEYDOWN into WM_CHAR, including for Backspace/Enter/Esc
        // -- see UiInput::textInput's own comment for why that is enough.
        if (wParam > 0 && wParam < 0x80) {
            g_workspacePendingText.push_back(static_cast<char>(wParam));
        }
        return 0;
    default:
        return DefWindowProc(hwnd, message, wParam, lParam);
    }
}
} // namespace

// Hosts every mode inside one window, one GL context and one gl33::Ui,
// switching between them by calling shutdown()/init() rather than
// destroying and recreating the window -- see Mode.hpp's own header
// comment for why (this replaces the earlier launcher_mode, which reopened
// a window per mode and made the WM_QUIT-drain bugfix necessary; a single
// persistent window removes that whole hazard class structurally instead
// of patching around it).
class Workspace {
public:
    explicit Workspace(std::string meshPath) {
        modes_.push_back(std::make_unique<MeshMode>(std::move(meshPath)));
        modes_.push_back(std::make_unique<HeatmapMode>());
        modes_.push_back(std::make_unique<CavityMode>());
        modes_.push_back(std::make_unique<Cavity3DMode>());
        modes_.push_back(std::make_unique<IsosurfaceMode>());
        modes_.push_back(std::make_unique<TurbulenceMode>());
        modes_.push_back(std::make_unique<SimMode>());
        modes_.push_back(std::make_unique<Sim3DMode>());
    }

    int run(int startIndex) {
        std::printf("Aether Workspace - navegue entre os modos pela barra lateral\n\n");
        std::fflush(stdout);

        g_workspaceWindow = gl33::createSimpleWindow(L"AetherWorkspace", L"Aether Workspace", g_workspaceWidth,
                                                       g_workspaceHeight, workspaceWndProc);
        if (!g_workspaceWindow) {
            std::fprintf(stderr, "erro ao criar janela\n");
            return 1;
        }
        HDC hdc = nullptr;
        HGLRC context = gl33::createGl33Context(g_workspaceWindow, hdc);
        if (!context) {
            return 1;
        }

        gl33::Ui ui;
        if (!ui.initialize()) {
            std::fprintf(stderr, "erro ao inicializar a UI\n");
            return 1;
        }

        currentIndex_ = std::clamp(startIndex, 0, static_cast<int>(modes_.size()) - 1);
        if (!modes_[currentIndex_]->init()) {
            std::fprintf(stderr, "erro ao inicializar o modo '%s'\n", modes_[currentIndex_]->label());
            return 1;
        }

        ShowWindow(g_workspaceWindow, SW_SHOW);
        UpdateWindow(g_workspaceWindow);

        // Same phantom-click hazard the old launcher_mode found and fixed
        // (a WM_LBUTTONDOWN already queued for whatever previously occupied
        // this window's screen position can land on the very first
        // rendered frame): a single warm-up here, since this window is now
        // created exactly once for the whole app lifetime rather than once
        // per mode switch.
        constexpr int kInputWarmupFrames = 15;
        int framesSinceOpen = 0;

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

            ++framesSinceOpen;
            const bool warmedUp = framesSinceOpen > kInputWarmupFrames;
            g_workspaceInput.mousePressed = warmedUp && g_workspacePendingPress;
            g_workspaceInput.mouseReleased = warmedUp && g_workspacePendingRelease;
            g_workspaceInput.textInput = g_workspacePendingText;
            const double wheelDelta = warmedUp ? g_workspacePendingWheelDelta : 0.0;
            g_workspacePendingPress = false;
            g_workspacePendingRelease = false;
            g_workspacePendingText.clear();
            g_workspacePendingWheelDelta = 0.0;

            ui.begin(g_workspaceWidth, g_workspaceHeight, g_workspaceInput);

            ui.beginPanel(kSidebarMargin, kSidebarMargin, kModeSwitcherWidth, "AETHER");
            for (std::size_t i = 0; i < modes_.size(); ++i) {
                std::string caption = (static_cast<int>(i) == currentIndex_) ? "[x] " : "[ ] ";
                caption += modes_[i]->label();
                if (ui.button(caption) && static_cast<int>(i) != currentIndex_) {
                    modes_[currentIndex_]->shutdown();
                    currentIndex_ = static_cast<int>(i);
                    if (!modes_[currentIndex_]->init()) {
                        std::fprintf(stderr, "erro ao inicializar o modo '%s'\n",
                                     modes_[currentIndex_]->label());
                        quit = true;
                    }
                }
            }
            ui.separator();
            ui.label("MALHA sem argumento:");
            ui.label("icosaedro em memoria. Para");
            ui.label("um STL: mesh <arquivo> pelo");
            ui.label("terminal.");
            ui.endPanel();

            if (!quit) {
                Mode& mode = *modes_[currentIndex_];
                mode.handleInput(g_workspaceInput, ui.wantsMouse(), wheelDelta);
                mode.update();
                mode.renderScene(g_workspaceWidth, g_workspaceHeight);
                mode.renderPanel(ui);
            }

            ui.end();

            SwapBuffers(hdc);
            Sleep(1);
        }

        modes_[currentIndex_]->shutdown();

        ui.shutdown();
        wglMakeCurrent(nullptr, nullptr);
        wglDeleteContext(context);
        ReleaseDC(g_workspaceWindow, hdc);
        DestroyWindow(g_workspaceWindow);
        g_workspaceWindow = nullptr;
        return 0;
    }

private:
    std::vector<std::unique_ptr<Mode>> modes_;
    int currentIndex_ = 0;
};

} // namespace aether_app_modes


int main(int argc, char** argv) {
    // Mode index order matches aether_app_modes::Workspace's constructor:
    // 0=mesh, 1=heatmap, 2=cavity, 3=cavity3d, 4=isosurface, 5=turbulence,
    // 6=sim, 7=sim3d. argv[1] (or its absence) only picks which one starts
    // selected -- the sidebar is always visible, so every mode is one click
    // away regardless of the starting choice.
    int startIndex = 0;
    std::string meshPath;

    if (argc >= 2) {
        const char* mode = argv[1];
        if (std::strcmp(mode, "hub") == 0 || std::strcmp(mode, "launcher") == 0) {
            startIndex = 0;
        } else if (std::strcmp(mode, "mesh") == 0) {
            startIndex = 0;
            if (argc >= 3) {
                meshPath = argv[2];
            }
        } else if (std::strcmp(mode, "heatmap") == 0) {
            startIndex = 1;
        } else if (std::strcmp(mode, "cavity") == 0) {
            startIndex = 2;
        } else if (std::strcmp(mode, "cavity3d") == 0) {
            startIndex = 3;
        } else if (std::strcmp(mode, "isosurface") == 0) {
            startIndex = 4;
        } else if (std::strcmp(mode, "turbulence") == 0) {
            startIndex = 5;
        } else if (std::strcmp(mode, "sim") == 0) {
            startIndex = 6;
        } else if (std::strcmp(mode, "sim3d") == 0) {
            startIndex = 7;
        } else {
            std::fprintf(stderr, "modo desconhecido: %s\n", mode);
            std::fprintf(stderr, "modos disponiveis: hub, mesh, heatmap, cavity, cavity3d, isosurface, "
                                  "turbulence, sim, sim3d\n");
            std::fprintf(stderr, "sem argumento nenhum abre o workspace (sidebar com todos os modos)\n");
            return 1;
        }
    }

    aether_app_modes::Workspace workspace(std::move(meshPath));
    return workspace.run(startIndex);
}
