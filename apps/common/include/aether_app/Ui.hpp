#pragma once

// Module 9 (CAD-like UI), foundation layer: a minimal immediate-mode GUI
// drawn with the same OpenGL 3.3 core-profile pipeline every viewer here
// already uses.
//
// **Why from scratch rather than Dear ImGui**, which is the obvious answer
// for exactly this problem: this project has a standing, load-bearing
// decision to have no external runtime dependencies -- apps/common exists
// precisely because the WGL/GL-3.3 bootstrap was hand-written instead of
// pulling in GLFW/GLAD (see Gl33.hpp). Adding an external UI library would
// reverse that decision, so it is not taken unilaterally here. The cost is
// real and worth naming: what follows is a few hundred lines that ImGui
// would have provided for free, with far fewer widgets. If the dependency
// stance is ever revisited, this is the layer ImGui would replace.
//
// **Design: immediate mode.** No retained widget tree, no callbacks, no
// invalidation -- each frame the caller calls panel()/label()/button()/
// slider()/checkbox() in order and gets the interaction result back
// immediately. Suits a tool UI whose state lives in the simulation
// objects, not in the widgets.
//
// **Rendering: one texture, one shader, one draw call per frame.** All
// geometry (panel backgrounds, button faces, slider tracks, and every text
// glyph) accumulates into a single vertex buffer of screen-space
// position + atlas UV + RGBA, flushed once in end(). Solid rectangles
// sample a deliberately all-white cell in the font atlas, so they need no
// separate shader or draw path -- which is why there is exactly one of
// each.
//
// **The font** is a 5x7 bitmap covering all 95 printable ASCII characters
// (32..126), embedded as source data so there is no font file to ship or
// load. Its glyph table was authored and then *verified glyph by glyph* by
// rendering the whole set as ASCII art before being converted to C++ (and
// converted programmatically, not transcribed by hand -- transcribing 665
// numbers is exactly the kind of step that introduces a wrong pixel nobody
// notices until a specific letter looks odd).

#include "aether_app/Gl33.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace aether::app {

// Per-frame input snapshot. The caller owns collecting this from its own
// WndProc -- Ui deliberately does not install any window hook.
struct UiInput {
    int mouseX = 0;
    int mouseY = 0;
    bool mouseDown = false;     // button currently held
    bool mousePressed = false;  // went down this frame
    bool mouseReleased = false; // came up this frame
    // Raw WM_CHAR characters collected since the previous begin() call, for
    // textField(). Backspace (0x08), Enter (0x0D) and Escape (0x1B) arrive
    // here too -- Windows generates WM_CHAR for control characters as well
    // as printable ones once TranslateMessage() has seen the keydown, so a
    // single message handler covers typing and the three control keys
    // textField() needs without a separate WM_KEYDOWN/VK_* path.
    std::string textInput;
};

struct UiColor {
    float r = 1.0f;
    float g = 1.0f;
    float b = 1.0f;
    float a = 1.0f;
};

class Ui {
public:
    // Builds the atlas texture, shader and buffers. Requires a current GL
    // 3.3 context (createGl33Context() already called). Returns false if
    // any GL object failed to build.
    bool initialize();
    void shutdown();

    // Starts a frame. screenWidth/Height are the framebuffer size in
    // pixels; all coordinates below are in those pixels, y down from the
    // top-left (the natural convention for a UI, and the one Win32 mouse
    // coordinates already arrive in).
    void begin(int screenWidth, int screenHeight, const UiInput& input);
    // Flushes every accumulated vertex in a single draw call.
    void end();

    // Opens a panel and makes it the layout target: subsequent widgets
    // stack downward inside it. Only one panel is open at a time.
    void beginPanel(int x, int y, int width, const char* title);
    void endPanel();

    // Widgets. Each returns whether the user interacted with it this
    // frame; slider/checkbox also write through their value pointer.
    void label(const std::string& text);
    void separator();
    bool button(const std::string& text);
    bool slider(const std::string& text, double* value, double minimum, double maximum);
    bool checkbox(const std::string& text, bool* value);
    // Click-to-focus numeric entry: shows "<label>: <value>" then an
    // editable box. While focused, digits/'.'/'-' type into a scratch
    // buffer, Backspace edits it, Enter commits (parses, clamps to
    // [minimum,maximum], writes *value, returns true that frame), Escape
    // discards, and a click anywhere outside the box also commits -- the
    // usual "click away" convention for a text field, not just a slider's
    // drag-while-held one. An unparseable buffer on commit (e.g. empty, or
    // just "-") leaves *value unchanged rather than writing garbage.
    bool textField(const std::string& text, double* value, double minimum, double maximum);

    // Direct drawing, for callers that want text or a rectangle outside
    // the panel layout (status lines, overlays).
    void drawText(int x, int y, const std::string& text, const UiColor& color);
    void drawRect(int x, int y, int width, int height, const UiColor& color);

    // Glyph cell metrics, so callers can size their own layouts.
    static constexpr int kGlyphWidth = 5;
    static constexpr int kGlyphHeight = 7;
    static constexpr int kGlyphAdvance = 6; // one blank column between glyphs
    static constexpr int kLineHeight = 14;

    int textWidth(const std::string& text) const {
        return static_cast<int>(text.size()) * kGlyphAdvance;
    }
    // True while the pointer is over the currently open panel -- lets the
    // host app suppress camera dragging when the click belongs to the UI.
    bool wantsMouse() const { return wantsMouse_; }

private:
    void pushQuad(float x, float y, float w, float h, float u0, float v0, float u1, float v1,
                   const UiColor& color);
    bool hit(int x, int y, int w, int h) const;
    // Parses editBuffer_, clamps into range and writes *value on success;
    // always clears editingField_. Shared by textField()'s three commit
    // paths (Enter, click-away, and -- trivially -- never leaving it
    // uncommitted) so the parse/clamp logic exists exactly once.
    bool commitEdit(double* value, double minimum, double maximum);

    GLuint program_ = 0;
    GLuint vao_ = 0;
    GLuint vbo_ = 0;
    GLuint atlas_ = 0;
    GLint screenUniform_ = -1;

    std::vector<float> vertices_;
    int screenWidth_ = 0;
    int screenHeight_ = 0;
    UiInput input_;
    bool wantsMouse_ = false;

    // Layout cursor for the open panel.
    int panelX_ = 0;
    int panelY_ = 0;
    int panelWidth_ = 0;
    int cursorY_ = 0;
    std::size_t panelBackgroundVertex_ = 0; // patched in endPanel() once the height is known
    bool panelOpen_ = false;
    // Identifies which slider is being dragged across frames; immediate
    // mode still needs this one piece of retained state, since a drag
    // spans frames by definition.
    int activeSlider_ = -1;
    int sliderCounter_ = 0;
    // Which textField() (by call-order index, matching activeSlider_'s own
    // convention) is currently being typed into, and its scratch buffer.
    // -1 means no field is focused. Deliberately NOT reset in begin() --
    // unlike activeSlider_, an edit-in-progress must survive frames where
    // the mouse isn't held at all.
    int editingField_ = -1;
    std::string editBuffer_;
    int fieldCounter_ = 0;
};

} // namespace aether::app
