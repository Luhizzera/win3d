#include "aether_app/Gl33.hpp"

#include <cstdio>

namespace aether::app {

PfnGlGenVertexArrays glGenVertexArrays = nullptr;
PfnGlBindVertexArray glBindVertexArray = nullptr;
PfnGlDeleteVertexArrays glDeleteVertexArrays = nullptr;
PfnGlGenBuffers glGenBuffers = nullptr;
PfnGlBindBuffer glBindBuffer = nullptr;
PfnGlBufferData glBufferData = nullptr;
PfnGlDeleteBuffers glDeleteBuffers = nullptr;
PfnGlVertexAttribPointer glVertexAttribPointer = nullptr;
PfnGlEnableVertexAttribArray glEnableVertexAttribArray = nullptr;
PfnGlCreateShader glCreateShader = nullptr;
PfnGlShaderSource glShaderSource = nullptr;
PfnGlCompileShader glCompileShader = nullptr;
PfnGlGetShaderiv glGetShaderiv = nullptr;
PfnGlGetShaderInfoLog glGetShaderInfoLog = nullptr;
PfnGlDeleteShader glDeleteShader = nullptr;
PfnGlCreateProgram glCreateProgram = nullptr;
PfnGlAttachShader glAttachShader = nullptr;
PfnGlLinkProgram glLinkProgram = nullptr;
PfnGlGetProgramiv glGetProgramiv = nullptr;
PfnGlGetProgramInfoLog glGetProgramInfoLog = nullptr;
PfnGlUseProgram glUseProgram = nullptr;
PfnGlDeleteProgram glDeleteProgram = nullptr;
PfnGlGetUniformLocation glGetUniformLocation = nullptr;
PfnGlUniformMatrix4fv glUniformMatrix4fv = nullptr;
PfnGlUniform3f glUniform3f = nullptr;
PfnGlUniform1f glUniform1f = nullptr;
PfnGlUniform2f glUniform2f = nullptr;

namespace {

constexpr int kWglContextMajorVersionArb = 0x2091;
constexpr int kWglContextMinorVersionArb = 0x2092;
constexpr int kWglContextProfileMaskArb = 0x9126;
constexpr int kWglContextCoreProfileBitArb = 0x00000001;

using PfnWglCreateContextAttribsArb = HGLRC(WINAPI*)(HDC, HGLRC, const int*);

template <typename Fn>
bool loadProc(const char* name, Fn& out) {
    out = reinterpret_cast<Fn>(wglGetProcAddress(name));
    return out != nullptr;
}

bool loadGl33Functions() {
    bool ok = true;
    ok &= loadProc("glGenVertexArrays", glGenVertexArrays);
    ok &= loadProc("glBindVertexArray", glBindVertexArray);
    ok &= loadProc("glDeleteVertexArrays", glDeleteVertexArrays);
    ok &= loadProc("glGenBuffers", glGenBuffers);
    ok &= loadProc("glBindBuffer", glBindBuffer);
    ok &= loadProc("glBufferData", glBufferData);
    ok &= loadProc("glDeleteBuffers", glDeleteBuffers);
    ok &= loadProc("glVertexAttribPointer", glVertexAttribPointer);
    ok &= loadProc("glEnableVertexAttribArray", glEnableVertexAttribArray);
    ok &= loadProc("glCreateShader", glCreateShader);
    ok &= loadProc("glShaderSource", glShaderSource);
    ok &= loadProc("glCompileShader", glCompileShader);
    ok &= loadProc("glGetShaderiv", glGetShaderiv);
    ok &= loadProc("glGetShaderInfoLog", glGetShaderInfoLog);
    ok &= loadProc("glDeleteShader", glDeleteShader);
    ok &= loadProc("glCreateProgram", glCreateProgram);
    ok &= loadProc("glAttachShader", glAttachShader);
    ok &= loadProc("glLinkProgram", glLinkProgram);
    ok &= loadProc("glGetProgramiv", glGetProgramiv);
    ok &= loadProc("glGetProgramInfoLog", glGetProgramInfoLog);
    ok &= loadProc("glUseProgram", glUseProgram);
    ok &= loadProc("glDeleteProgram", glDeleteProgram);
    ok &= loadProc("glGetUniformLocation", glGetUniformLocation);
    ok &= loadProc("glUniformMatrix4fv", glUniformMatrix4fv);
    ok &= loadProc("glUniform3f", glUniform3f);
    ok &= loadProc("glUniform1f", glUniform1f);
    ok &= loadProc("glUniform2f", glUniform2f);
    return ok;
}

} // namespace

PIXELFORMATDESCRIPTOR standardPixelFormatDescriptor() {
    PIXELFORMATDESCRIPTOR pfd{};
    pfd.nSize = sizeof(pfd);
    pfd.nVersion = 1;
    pfd.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER;
    pfd.iPixelType = PFD_TYPE_RGBA;
    pfd.cColorBits = 32;
    pfd.cDepthBits = 24;
    pfd.iLayerType = PFD_MAIN_PLANE;
    return pfd;
}

HWND createSimpleWindow(const wchar_t* className, const wchar_t* title, int width, int height,
                         WNDPROC wndProc) {
    WNDCLASSEXW wc{};
    wc.cbSize = sizeof(wc);
    wc.style = CS_OWNDC;
    wc.lpfnWndProc = wndProc;
    wc.hInstance = GetModuleHandleW(nullptr);
    wc.hCursor = LoadCursorW(nullptr, IDC_ARROW);
    wc.lpszClassName = className;
    RegisterClassExW(&wc);

    return CreateWindowExW(0, className, title, WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, width,
                            height, nullptr, nullptr, wc.hInstance, nullptr);
}

HGLRC createGl33Context(HWND hwnd, HDC& hOut) {
    // Step 1: throwaway legacy-context window, purely to fetch
    // wglCreateContextAttribsARB. See the header comment for why it must
    // use DefWindowProc, not the caller's real WndProc.
    HWND dummyHwnd = createSimpleWindow(L"AetherAppGl33DummyWindowClass", L"", 4, 4, DefWindowProc);
    HDC dummyHdc = GetDC(dummyHwnd);
    PIXELFORMATDESCRIPTOR dummyPfd = standardPixelFormatDescriptor();
    const int dummyPixelFormat = ChoosePixelFormat(dummyHdc, &dummyPfd);
    SetPixelFormat(dummyHdc, dummyPixelFormat, &dummyPfd);
    HGLRC dummyContext = wglCreateContext(dummyHdc);
    wglMakeCurrent(dummyHdc, dummyContext);

    PfnWglCreateContextAttribsArb wglCreateContextAttribsARB = nullptr;
    loadProc("wglCreateContextAttribsARB", wglCreateContextAttribsARB);

    wglMakeCurrent(nullptr, nullptr);
    wglDeleteContext(dummyContext);
    ReleaseDC(dummyHwnd, dummyHdc);
    DestroyWindow(dummyHwnd);

    if (!wglCreateContextAttribsARB) {
        std::fprintf(stderr, "wglCreateContextAttribsARB indisponivel - driver sem suporte a GL 3.3 core\n");
        return nullptr;
    }

    // Step 2: the real context, on the caller's real window.
    HDC hdc = GetDC(hwnd);
    PIXELFORMATDESCRIPTOR pfd = standardPixelFormatDescriptor();
    const int pixelFormat = ChoosePixelFormat(hdc, &pfd);
    SetPixelFormat(hdc, pixelFormat, &pfd);

    const int contextAttribs[] = {kWglContextMajorVersionArb,
                                   3,
                                   kWglContextMinorVersionArb,
                                   3,
                                   kWglContextProfileMaskArb,
                                   kWglContextCoreProfileBitArb,
                                   0};
    HGLRC context = wglCreateContextAttribsARB(hdc, nullptr, contextAttribs);
    if (!context) {
        std::fprintf(stderr, "erro ao criar contexto OpenGL 3.3 core\n");
        ReleaseDC(hwnd, hdc);
        return nullptr;
    }
    wglMakeCurrent(hdc, context);

    if (!loadGl33Functions()) {
        std::fprintf(stderr, "erro ao carregar funcoes OpenGL 3.3\n");
        wglMakeCurrent(nullptr, nullptr);
        wglDeleteContext(context);
        ReleaseDC(hwnd, hdc);
        return nullptr;
    }

    hOut = hdc;
    return context;
}

GLuint compileShader(GLenum type, const char* source) {
    const GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    GLint success = 0;
    glGetShaderiv(shader, kGlCompileStatus, &success);
    if (!success) {
        char log[1024];
        GLsizei logLength = 0;
        glGetShaderInfoLog(shader, sizeof(log), &logLength, log);
        std::fprintf(stderr, "erro ao compilar shader: %.*s\n", logLength, log);
    }
    return shader;
}

GLuint linkProgram(GLuint vertexShader, GLuint fragmentShader) {
    const GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);
    GLint success = 0;
    glGetProgramiv(program, kGlLinkStatus, &success);
    if (!success) {
        char log[1024];
        GLsizei logLength = 0;
        glGetProgramInfoLog(program, sizeof(log), &logLength, log);
        std::fprintf(stderr, "erro ao linkar shader program: %.*s\n", logLength, log);
    }
    return program;
}

} // namespace aether::app
