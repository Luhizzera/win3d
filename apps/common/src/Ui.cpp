#include "aether_app/Ui.hpp"

#include <algorithm>
#include <cstdio>

namespace aether::app {

namespace {

// GL constants <GL/gl.h> (OpenGL 1.1) does not define but the 3.3 core
// profile needs. GL_LUMINANCE and GL_ALPHA were removed in core profile, so
// the single-channel atlas uses GL_RED / GL_R8 instead.
constexpr GLenum kGlRed = 0x1903;
constexpr GLenum kGlR8 = 0x8229;
constexpr GLenum kGlClampToEdge = 0x812F;

// 5x7 bitmap font covering printable ASCII 32..126. Each row's low 5 bits
// are the pixels, the most significant of those five being the leftmost.
//
// **Generated, not hand-transcribed**, from a table that was verified glyph
// by glyph by rendering the whole set as ASCII art first (see Ui.hpp).
// Transcribing 665 numbers by hand is exactly the step that plants a wrong
// pixel nobody notices until one specific letter looks subtly off. If a
// glyph ever needs changing, change it in a form you can *see* and
// regenerate.
constexpr int kFirstGlyph = 32;
constexpr int kGlyphCount = 95;
constexpr unsigned char kFont[kGlyphCount][7] = {
    { 0,  0,  0,  0,  0,  0,  0}, // ' '
    { 4,  4,  4,  4,  4,  0,  4}, // '!'
    {10, 10, 10,  0,  0,  0,  0}, // '"'
    {10, 10, 31, 10, 31, 10, 10}, // '#'
    { 4, 15, 20, 14,  5, 30,  4}, // '$'
    {24, 25,  2,  4,  8, 19,  3}, // '%'
    { 8, 20, 20,  8, 21, 18, 13}, // '&'
    {12,  4,  8,  0,  0,  0,  0}, // "'"
    { 2,  4,  8,  8,  8,  4,  2}, // '('
    { 8,  4,  2,  2,  2,  4,  8}, // ')'
    { 0,  4, 21, 14, 21,  4,  0}, // '*'
    { 0,  4,  4, 31,  4,  4,  0}, // '+'
    { 0,  0,  0,  0,  0,  4,  8}, // ','
    { 0,  0,  0, 31,  0,  0,  0}, // '-'
    { 0,  0,  0,  0,  0, 12, 12}, // '.'
    { 0,  1,  2,  4,  8, 16,  0}, // '/'
    {14, 17, 19, 21, 25, 17, 14}, // '0'
    { 4, 12,  4,  4,  4,  4, 14}, // '1'
    {14, 17,  1,  2,  4,  8, 31}, // '2'
    {31,  2,  4,  2,  1, 17, 14}, // '3'
    { 2,  6, 10, 18, 31,  2,  2}, // '4'
    {31, 16, 30,  1,  1, 17, 14}, // '5'
    { 6,  8, 16, 30, 17, 17, 14}, // '6'
    {31,  1,  2,  4,  8,  8,  8}, // '7'
    {14, 17, 17, 14, 17, 17, 14}, // '8'
    {14, 17, 17, 15,  1,  2, 12}, // '9'
    { 0, 12, 12,  0, 12, 12,  0}, // ':'
    { 0, 12, 12,  0, 12,  4,  8}, // ';'
    { 2,  4,  8, 16,  8,  4,  2}, // '<'
    { 0,  0, 31,  0, 31,  0,  0}, // '='
    { 8,  4,  2,  1,  2,  4,  8}, // '>'
    {14, 17,  1,  2,  4,  0,  4}, // '?'
    {14, 17,  1, 13, 21, 21, 14}, // '@'
    { 4, 10, 17, 17, 31, 17, 17}, // 'A'
    {30, 17, 17, 30, 17, 17, 30}, // 'B'
    {14, 17, 16, 16, 16, 17, 14}, // 'C'
    {28, 18, 17, 17, 17, 18, 28}, // 'D'
    {31, 16, 16, 30, 16, 16, 31}, // 'E'
    {31, 16, 16, 30, 16, 16, 16}, // 'F'
    {14, 17, 16, 23, 17, 17, 15}, // 'G'
    {17, 17, 17, 31, 17, 17, 17}, // 'H'
    {14,  4,  4,  4,  4,  4, 14}, // 'I'
    { 7,  2,  2,  2,  2, 18, 12}, // 'J'
    {17, 18, 20, 24, 20, 18, 17}, // 'K'
    {16, 16, 16, 16, 16, 16, 31}, // 'L'
    {17, 27, 21, 21, 17, 17, 17}, // 'M'
    {17, 17, 25, 21, 19, 17, 17}, // 'N'
    {14, 17, 17, 17, 17, 17, 14}, // 'O'
    {30, 17, 17, 30, 16, 16, 16}, // 'P'
    {14, 17, 17, 17, 21, 18, 13}, // 'Q'
    {30, 17, 17, 30, 20, 18, 17}, // 'R'
    {15, 16, 16, 14,  1,  1, 30}, // 'S'
    {31,  4,  4,  4,  4,  4,  4}, // 'T'
    {17, 17, 17, 17, 17, 17, 14}, // 'U'
    {17, 17, 17, 17, 17, 10,  4}, // 'V'
    {17, 17, 17, 21, 21, 21, 10}, // 'W'
    {17, 17, 10,  4, 10, 17, 17}, // 'X'
    {17, 17, 10,  4,  4,  4,  4}, // 'Y'
    {31,  1,  2,  4,  8, 16, 31}, // 'Z'
    {14,  8,  8,  8,  8,  8, 14}, // '['
    { 0, 16,  8,  4,  2,  1,  0}, // '\\'
    {14,  2,  2,  2,  2,  2, 14}, // ']'
    { 4, 10, 17,  0,  0,  0,  0}, // '^'
    { 0,  0,  0,  0,  0,  0, 31}, // '_'
    { 8,  4,  2,  0,  0,  0,  0}, // '`'
    { 0,  0, 14,  1, 15, 17, 15}, // 'a'
    {16, 16, 30, 17, 17, 17, 30}, // 'b'
    { 0,  0, 14, 17, 16, 17, 14}, // 'c'
    { 1,  1, 15, 17, 17, 17, 15}, // 'd'
    { 0,  0, 14, 17, 31, 16, 14}, // 'e'
    { 6,  9,  8, 28,  8,  8,  8}, // 'f'
    { 0,  0, 15, 17, 15,  1, 14}, // 'g'
    {16, 16, 30, 17, 17, 17, 17}, // 'h'
    { 4,  0, 12,  4,  4,  4, 14}, // 'i'
    { 2,  0,  6,  2,  2, 18, 12}, // 'j'
    {16, 16, 18, 20, 24, 20, 18}, // 'k'
    {12,  4,  4,  4,  4,  4, 14}, // 'l'
    { 0,  0, 26, 21, 21, 21, 21}, // 'm'
    { 0,  0, 30, 17, 17, 17, 17}, // 'n'
    { 0,  0, 14, 17, 17, 17, 14}, // 'o'
    { 0,  0, 30, 17, 30, 16, 16}, // 'p'
    { 0,  0, 15, 17, 15,  1,  1}, // 'q'
    { 0,  0, 22, 25, 16, 16, 16}, // 'r'
    { 0,  0, 15, 16, 14,  1, 30}, // 's'
    { 8,  8, 28,  8,  8,  9,  6}, // 't'
    { 0,  0, 17, 17, 17, 19, 13}, // 'u'
    { 0,  0, 17, 17, 17, 10,  4}, // 'v'
    { 0,  0, 17, 17, 21, 21, 10}, // 'w'
    { 0,  0, 17, 10,  4, 10, 17}, // 'x'
    { 0,  0, 17, 17, 15,  1, 14}, // 'y'
    { 0,  0, 31,  2,  4,  8, 31}, // 'z'
    { 2,  4,  4,  8,  4,  4,  2}, // '{'
    { 4,  4,  4,  4,  4,  4,  4}, // '|'
    { 8,  4,  4,  2,  4,  4,  8}, // '}'
    { 0,  0,  8, 21,  2,  0,  0}, // '~'
};

// Atlas layout: 16 glyph cells per row, 6 rows = 96 cells for 95 glyphs.
// The one spare cell (index 95, last) is filled solid white and is what
// every non-text rectangle samples -- which is why panels, buttons, slider
// tracks and text all go through one shader and one draw call.
constexpr int kCellsPerRow = 16;
constexpr int kCellRows = 6;
constexpr int kCellWidth = 6;  // 5 glyph columns + 1 spacer
constexpr int kCellHeight = 8; // 7 glyph rows + 1 spacer
constexpr int kAtlasWidth = kCellsPerRow * kCellWidth;   // 96
constexpr int kAtlasHeight = kCellRows * kCellHeight;    // 48
constexpr int kSolidCell = 95;

const char* kVertexShader = R"(#version 330 core
layout(location = 0) in vec2 aPos;
layout(location = 1) in vec2 aUv;
layout(location = 2) in vec4 aColor;
uniform vec2 uScreen;
out vec2 vUv;
out vec4 vColor;
void main() {
    // Screen pixels, y down from top-left, into clip space.
    vec2 ndc = vec2(aPos.x / uScreen.x * 2.0 - 1.0, 1.0 - aPos.y / uScreen.y * 2.0);
    gl_Position = vec4(ndc, 0.0, 1.0);
    vUv = aUv;
    vColor = aColor;
}
)";

const char* kFragmentShader = R"(#version 330 core
in vec2 vUv;
in vec4 vColor;
uniform sampler2D uAtlas;
out vec4 FragColor;
void main() {
    // The atlas is a coverage mask: solid shapes sample the all-white cell
    // and come out at full alpha, glyph pixels mask the text colour.
    float coverage = texture(uAtlas, vUv).r;
    FragColor = vec4(vColor.rgb, vColor.a * coverage);
}
)";

// Palette -- a muted CAD-ish dark theme.
const UiColor kPanelBackground{0.12f, 0.13f, 0.15f, 0.92f};
const UiColor kPanelTitle{0.55f, 0.78f, 1.00f, 1.0f};
const UiColor kText{0.86f, 0.88f, 0.90f, 1.0f};
const UiColor kTextDim{0.55f, 0.58f, 0.62f, 1.0f};
const UiColor kButtonFace{0.22f, 0.24f, 0.28f, 1.0f};
const UiColor kButtonHover{0.30f, 0.34f, 0.40f, 1.0f};
const UiColor kAccent{0.30f, 0.62f, 0.95f, 1.0f};
const UiColor kTrack{0.18f, 0.19f, 0.22f, 1.0f};

constexpr int kPadding = 8;
constexpr int kWidgetHeight = 18;
constexpr int kTitleHeight = 22;

} // namespace

bool Ui::initialize() {
    // --- atlas ---
    std::vector<unsigned char> pixels(static_cast<std::size_t>(kAtlasWidth) * kAtlasHeight, 0);
    for (int g = 0; g < kGlyphCount; ++g) {
        const int cellX = (g % kCellsPerRow) * kCellWidth;
        const int cellY = (g / kCellsPerRow) * kCellHeight;
        for (int row = 0; row < 7; ++row) {
            const unsigned char bits = kFont[g][row];
            for (int col = 0; col < 5; ++col) {
                if ((bits & (1u << (4 - col))) != 0u) {
                    pixels[static_cast<std::size_t>(cellY + row) * kAtlasWidth + (cellX + col)] = 255;
                }
            }
        }
    }
    {
        const int cellX = (kSolidCell % kCellsPerRow) * kCellWidth;
        const int cellY = (kSolidCell / kCellsPerRow) * kCellHeight;
        for (int row = 0; row < kCellHeight; ++row) {
            for (int col = 0; col < kCellWidth; ++col) {
                pixels[static_cast<std::size_t>(cellY + row) * kAtlasWidth + (cellX + col)] = 255;
            }
        }
    }

    glGenTextures(1, &atlas_);
    glBindTexture(GL_TEXTURE_2D, atlas_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // rows are 96 bytes but not 4-aligned in general
    glTexImage2D(GL_TEXTURE_2D, 0, static_cast<GLint>(kGlR8), kAtlasWidth, kAtlasHeight, 0, kGlRed,
                 GL_UNSIGNED_BYTE, pixels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, static_cast<GLint>(kGlClampToEdge));
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, static_cast<GLint>(kGlClampToEdge));

    // --- shader ---
    const GLuint vertexShader = compileShader(kGlVertexShader, kVertexShader);
    const GLuint fragmentShader = compileShader(kGlFragmentShader, kFragmentShader);
    program_ = linkProgram(vertexShader, fragmentShader);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    if (program_ == 0) {
        return false;
    }
    screenUniform_ = glGetUniformLocation(program_, "uScreen");

    // --- buffers ---
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);
    glGenBuffers(1, &vbo_);
    glBindBuffer(kGlArrayBuffer, vbo_);
    constexpr GLsizei stride = 8 * sizeof(float);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(2 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void*>(4 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glBindVertexArray(0);
    return true;
}

void Ui::shutdown() {
    if (vbo_ != 0) {
        glDeleteBuffers(1, &vbo_);
        vbo_ = 0;
    }
    if (vao_ != 0) {
        glDeleteVertexArrays(1, &vao_);
        vao_ = 0;
    }
    if (program_ != 0) {
        glDeleteProgram(program_);
        program_ = 0;
    }
    if (atlas_ != 0) {
        glDeleteTextures(1, &atlas_);
        atlas_ = 0;
    }
}

void Ui::begin(int screenWidth, int screenHeight, const UiInput& input) {
    screenWidth_ = screenWidth;
    screenHeight_ = screenHeight;
    input_ = input;
    vertices_.clear();
    wantsMouse_ = false;
    sliderCounter_ = 0;
    if (!input_.mouseDown) {
        activeSlider_ = -1; // a drag cannot outlive the button being held
    }
}

void Ui::pushQuad(float x, float y, float w, float h, float u0, float v0, float u1, float v1,
                   const UiColor& c) {
    const float quad[6][4] = {{x, y, u0, v0},         {x + w, y, u1, v0},     {x + w, y + h, u1, v1},
                              {x, y, u0, v0},         {x + w, y + h, u1, v1}, {x, y + h, u0, v1}};
    for (const auto& v : quad) {
        vertices_.insert(vertices_.end(), {v[0], v[1], v[2], v[3], c.r, c.g, c.b, c.a});
    }
}

void Ui::drawRect(int x, int y, int width, int height, const UiColor& color) {
    // Sample the middle of the solid cell, well away from its edges, so
    // nearest-neighbour filtering can never pick up a neighbouring glyph.
    const float cx = (static_cast<float>((kSolidCell % kCellsPerRow) * kCellWidth) + 2.5f) / kAtlasWidth;
    const float cy = (static_cast<float>((kSolidCell / kCellsPerRow) * kCellHeight) + 3.5f) / kAtlasHeight;
    pushQuad(static_cast<float>(x), static_cast<float>(y), static_cast<float>(width),
             static_cast<float>(height), cx, cy, cx, cy, color);
}

void Ui::drawText(int x, int y, const std::string& text, const UiColor& color) {
    int penX = x;
    for (const char rawChar : text) {
        const int code = static_cast<unsigned char>(rawChar);
        if (code >= kFirstGlyph && code < kFirstGlyph + kGlyphCount) {
            const int g = code - kFirstGlyph;
            const float u0 = static_cast<float>((g % kCellsPerRow) * kCellWidth) / kAtlasWidth;
            const float v0 = static_cast<float>((g / kCellsPerRow) * kCellHeight) / kAtlasHeight;
            const float u1 = u0 + static_cast<float>(kGlyphWidth) / kAtlasWidth;
            const float v1 = v0 + static_cast<float>(kGlyphHeight) / kAtlasHeight;
            pushQuad(static_cast<float>(penX), static_cast<float>(y), static_cast<float>(kGlyphWidth),
                     static_cast<float>(kGlyphHeight), u0, v0, u1, v1, color);
        }
        penX += kGlyphAdvance;
    }
}

bool Ui::hit(int x, int y, int w, int h) const {
    return input_.mouseX >= x && input_.mouseX < x + w && input_.mouseY >= y && input_.mouseY < y + h;
}

void Ui::beginPanel(int x, int y, int width, const char* title) {
    panelX_ = x;
    panelY_ = y;
    panelWidth_ = width;
    panelOpen_ = true;
    // Reserve the background quad now and patch its height in endPanel(),
    // once the widgets have decided how tall the panel actually is.
    panelBackgroundVertex_ = vertices_.size();
    drawRect(x, y, width, kTitleHeight, kPanelBackground);
    drawText(x + kPadding, y + (kTitleHeight - kGlyphHeight) / 2, title, kPanelTitle);
    cursorY_ = y + kTitleHeight + kPadding;
}

void Ui::endPanel() {
    if (!panelOpen_) {
        return;
    }
    const int height = cursorY_ - panelY_ + kPadding - kWidgetHeight + kWidgetHeight;
    // Patch the reserved background quad's y-extent (vertex layout is
    // 6 vertices x 8 floats; y is component 1 of each).
    const float top = static_cast<float>(panelY_);
    const float bottom = static_cast<float>(panelY_ + height);
    const float ys[6] = {top, top, bottom, top, bottom, bottom};
    for (int i = 0; i < 6; ++i) {
        vertices_[panelBackgroundVertex_ + static_cast<std::size_t>(i) * 8 + 1] = ys[i];
    }
    // The title text was emitted after the background, so it still draws on
    // top; re-emit nothing here.
    wantsMouse_ = wantsMouse_ || hit(panelX_, panelY_, panelWidth_, height);
    panelOpen_ = false;
}

void Ui::label(const std::string& text) {
    drawText(panelX_ + kPadding, cursorY_ + (kWidgetHeight - kGlyphHeight) / 2, text, kText);
    cursorY_ += kLineHeight;
}

void Ui::separator() {
    cursorY_ += 4;
    drawRect(panelX_ + kPadding, cursorY_, panelWidth_ - 2 * kPadding, 1, kTextDim);
    cursorY_ += 8;
}

bool Ui::button(const std::string& text) {
    const int x = panelX_ + kPadding;
    const int y = cursorY_;
    const int w = panelWidth_ - 2 * kPadding;
    const bool hovered = hit(x, y, w, kWidgetHeight);
    drawRect(x, y, w, kWidgetHeight, hovered ? kButtonHover : kButtonFace);
    drawText(x + (w - textWidth(text)) / 2, y + (kWidgetHeight - kGlyphHeight) / 2, text, kText);
    cursorY_ += kWidgetHeight + 4;
    return hovered && input_.mousePressed;
}

bool Ui::checkbox(const std::string& text, bool* value) {
    const int x = panelX_ + kPadding;
    const int y = cursorY_;
    const int box = 12;
    const bool hovered = hit(x, y, panelWidth_ - 2 * kPadding, kWidgetHeight);
    const int boxY = y + (kWidgetHeight - box) / 2;
    drawRect(x, boxY, box, box, hovered ? kButtonHover : kButtonFace);
    if (*value) {
        drawRect(x + 3, boxY + 3, box - 6, box - 6, kAccent);
    }
    drawText(x + box + 6, y + (kWidgetHeight - kGlyphHeight) / 2, text, kText);
    cursorY_ += kWidgetHeight + 4;
    const bool clicked = hovered && input_.mousePressed;
    if (clicked) {
        *value = !*value;
    }
    return clicked;
}

bool Ui::slider(const std::string& text, double* value, double minimum, double maximum) {
    const int id = sliderCounter_++;
    const int x = panelX_ + kPadding;
    cursorY_ += 3; // a slider is two stacked rows, so it needs its own breathing room
    const int y = cursorY_;
    const int w = panelWidth_ - 2 * kPadding;

    char buffer[128];
    std::snprintf(buffer, sizeof(buffer), "%s: %.4g", text.c_str(), *value);
    drawText(x, y, buffer, kText);

    const int trackY = y + kGlyphHeight + 3;
    const int trackH = 8;
    drawRect(x, trackY, w, trackH, kTrack);

    const bool hovered = hit(x, trackY - 3, w, trackH + 6);
    if (hovered && input_.mousePressed) {
        activeSlider_ = id;
    }
    bool changed = false;
    if (activeSlider_ == id && input_.mouseDown) {
        const double t = std::clamp(static_cast<double>(input_.mouseX - x) / std::max(w, 1), 0.0, 1.0);
        const double updated = minimum + t * (maximum - minimum);
        if (updated != *value) {
            *value = updated;
            changed = true;
        }
    }

    const double span = (maximum - minimum) != 0.0 ? (maximum - minimum) : 1.0;
    const double fraction = std::clamp((*value - minimum) / span, 0.0, 1.0);
    drawRect(x, trackY, static_cast<int>(fraction * w), trackH, kAccent);
    const int knobX = x + static_cast<int>(fraction * w) - 2;
    drawRect(knobX, trackY - 2, 5, trackH + 4, kText);

    cursorY_ += kGlyphHeight + trackH + 10;
    return changed;
}

void Ui::end() {
    if (vertices_.empty()) {
        return;
    }
    // UI draws last, as a flat overlay: no depth, straight alpha blending.
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glUseProgram(program_);
    glUniform2f(screenUniform_, static_cast<GLfloat>(screenWidth_), static_cast<GLfloat>(screenHeight_));
    glBindTexture(GL_TEXTURE_2D, atlas_);
    glBindVertexArray(vao_);
    glBindBuffer(kGlArrayBuffer, vbo_);
    glBufferData(kGlArrayBuffer, static_cast<GLsizeiptr>(vertices_.size() * sizeof(float)), vertices_.data(),
                 kGlStaticDraw);
    glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(vertices_.size() / 8));
    glBindVertexArray(0);

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}

} // namespace aether::app
