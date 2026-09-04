#pragma once

#include "aether_app/Ui.hpp"

// The interface every viewer mode implements so `Workspace` (see main.cpp)
// can host all of them inside one window, one GL context and one
// `gl33::Ui` sidebar, switching between them with no window teardown --
// requested directly by the user after trying the window-per-mode
// `launcher_mode` hub and finding it unsatisfying (the hub disappeared
// while a mode was open, since each mode still blocked the thread running
// its own message loop until its own window closed).
//
// Every method below existed already in each mode's own `run()` before
// this interface was introduced -- this only draws the boundary between
// "runs once when a mode becomes active" (init), "runs once before a mode
// stops being active" (shutdown) and "runs every frame while active"
// (handleInput/update/renderScene/renderPanel), the same four phases every
// mode's `run()` already had inline, just not previously reusable across a
// mode switch that keeps the window and GL context alive.
//
// **`renderScene` draws full-window, `renderPanel` draws last.** Proven
// already by sim3d_mode before this interface existed: a 3D/2D scene
// rendered with the mode's own shader/VAO covering the whole framebuffer,
// then the Ui's own panel drawn on top as a flat overlay (Ui disables
// depth testing for its own draw), never a viewport carved out to make
// room for a sidebar. `Workspace::run()` calls renderScene() then
// renderPanel() inside the same ui.begin()/ui.end() pair, so a mode's
// panel content composes with the Workspace's own mode-switcher panel
// drawn just above it -- one `ui.end()` flush covers both.
namespace aether_app_modes {

class Mode {
public:
    virtual ~Mode() = default;

    // Button label in the Workspace's mode-switcher panel.
    virtual const char* label() const = 0;

    // Called once when this mode becomes the active one: compile shaders,
    // build VAO/VBO, run whatever one-time solve or bake every mode
    // already did before opening its own window. Returns false on failure
    // (mirrors every mode's existing `run()` returning 1 on setup failure).
    virtual bool init() = 0;

    // Called once when this mode stops being the active one (another mode
    // was chosen, or the Workspace is closing): delete this mode's own GL
    // objects. The window and GL context themselves are the Workspace's,
    // not this mode's, to tear down.
    virtual void shutdown() = 0;

    // `input` is this frame's snapshot (same aether_app::UiInput every
    // mode already builds from its own WndProc's queued press/release/text
    // -- see Ui.hpp's own comment on why Ui itself installs no window
    // hook). `uiWantsMouse` is `ui.wantsMouse()` from the frame's already-
    // drawn switcher panel, so a mode's own camera-drag logic can suppress
    // itself exactly the way sim3d_mode already does
    // (`if (input.mousePressed && !uiWantsMouse)`) without needing to know
    // panel geometry itself. `wheelDelta` is +1/-1 the one frame a
    // WM_MOUSEWHEEL arrived, 0 otherwise.
    virtual void handleInput(const aether::app::UiInput& input, bool uiWantsMouse, double wheelDelta) = 0;

    // Advances whatever this mode simulates, if anything -- a no-op for
    // every mode that only ever solved once up front (six of the eight
    // measured before this interface was written; see the migration's own
    // plan notes). Must never block: sim3d_mode's own worker-thread/
    // try_to_lock snapshot pattern is the template for a mode whose solve
    // is too slow to run inline per frame.
    virtual void update() = 0;

    // Draws this mode's own content across the full `width` x `height`
    // framebuffer, using its own shader/VAO/program. Called before the
    // Workspace's Ui overlay is flushed, never after.
    virtual void renderScene(int width, int height) = 0;

    // Draws this mode's own panel content (parameters, diagnostics -- can
    // be nothing at all for a mode with none) via `ui.beginPanel()`/
    // widgets/`ui.endPanel()`. Called between the Workspace's own
    // `ui.beginPanel()`/`ui.endPanel()` (the mode switcher) and `ui.end()`
    // -- never call `ui.begin()`/`ui.end()` here, the Workspace owns that
    // pairing for the whole frame.
    virtual void renderPanel(aether::app::Ui& ui) = 0;
};

} // namespace aether_app_modes
