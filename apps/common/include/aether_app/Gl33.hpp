#pragma once

// Shared OpenGL 3.3 core-profile + WGL bootstrap toolkit for every apps/*
// GUI viewer, extracted out of apps/viewer once it proved the pipeline
// works (see that app's own comment for the full rationale: no GLFW/GLAD/
// Vulkan, hand-loaded function pointers via wglGetProcAddress, the
// "throwaway legacy context" trick to fetch wglCreateContextAttribsARB).
// Centralizing this here, rather than copy-pasting apps/viewer's ~150
// lines of boilerplate into every other viewer, both avoids the
// duplication and closes off a real bug class already hit once (the
// WM_QUIT-via-shared-WndProc bug in the throwaway window -- see
// createGl33Context()'s own comment).

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#define UNICODE
#define _UNICODE

#include <windows.h>

#include <GL/gl.h>

#include <cstddef>

namespace aether::app {

using GLchar = char;
using GLsizeiptr = ptrdiff_t;

// GL 3.3 core constants <GL/gl.h> does not define (it only covers OpenGL
// 1.1).
constexpr GLenum kGlArrayBuffer = 0x8892;
constexpr GLenum kGlStaticDraw = 0x88E4;
constexpr GLenum kGlFragmentShader = 0x8B30;
constexpr GLenum kGlVertexShader = 0x8B31;
constexpr GLenum kGlCompileStatus = 0x8B81;
constexpr GLenum kGlLinkStatus = 0x8B82;

using PfnGlGenVertexArrays = void(WINAPI*)(GLsizei, GLuint*);
using PfnGlBindVertexArray = void(WINAPI*)(GLuint);
using PfnGlDeleteVertexArrays = void(WINAPI*)(GLsizei, const GLuint*);
using PfnGlGenBuffers = void(WINAPI*)(GLsizei, GLuint*);
using PfnGlBindBuffer = void(WINAPI*)(GLenum, GLuint);
using PfnGlBufferData = void(WINAPI*)(GLenum, GLsizeiptr, const void*, GLenum);
using PfnGlDeleteBuffers = void(WINAPI*)(GLsizei, const GLuint*);
using PfnGlVertexAttribPointer = void(WINAPI*)(GLuint, GLint, GLenum, GLboolean, GLsizei, const void*);
using PfnGlEnableVertexAttribArray = void(WINAPI*)(GLuint);
using PfnGlCreateShader = GLuint(WINAPI*)(GLenum);
using PfnGlShaderSource = void(WINAPI*)(GLuint, GLsizei, const GLchar* const*, const GLint*);
using PfnGlCompileShader = void(WINAPI*)(GLuint);
using PfnGlGetShaderiv = void(WINAPI*)(GLuint, GLenum, GLint*);
using PfnGlGetShaderInfoLog = void(WINAPI*)(GLuint, GLsizei, GLsizei*, GLchar*);
using PfnGlDeleteShader = void(WINAPI*)(GLuint);
using PfnGlCreateProgram = GLuint(WINAPI*)(void);
using PfnGlAttachShader = void(WINAPI*)(GLuint, GLuint);
using PfnGlLinkProgram = void(WINAPI*)(GLuint);
using PfnGlGetProgramiv = void(WINAPI*)(GLuint, GLenum, GLint*);
using PfnGlGetProgramInfoLog = void(WINAPI*)(GLuint, GLsizei, GLsizei*, GLchar*);
using PfnGlUseProgram = void(WINAPI*)(GLuint);
using PfnGlDeleteProgram = void(WINAPI*)(GLuint);
using PfnGlGetUniformLocation = GLint(WINAPI*)(GLuint, const GLchar*);
using PfnGlUniformMatrix4fv = void(WINAPI*)(GLint, GLsizei, GLboolean, const GLfloat*);
using PfnGlUniform3f = void(WINAPI*)(GLint, GLfloat, GLfloat, GLfloat);
using PfnGlUniform1f = void(WINAPI*)(GLint, GLfloat);

// Populated by createGl33Context() on success; every apps/* viewer shares
// these same pointers (one process, one GL context at a time -- true for
// every viewer in this project, each is its own single-window executable).
extern PfnGlGenVertexArrays glGenVertexArrays;
extern PfnGlBindVertexArray glBindVertexArray;
extern PfnGlDeleteVertexArrays glDeleteVertexArrays;
extern PfnGlGenBuffers glGenBuffers;
extern PfnGlBindBuffer glBindBuffer;
extern PfnGlBufferData glBufferData;
extern PfnGlDeleteBuffers glDeleteBuffers;
extern PfnGlVertexAttribPointer glVertexAttribPointer;
extern PfnGlEnableVertexAttribArray glEnableVertexAttribArray;
extern PfnGlCreateShader glCreateShader;
extern PfnGlShaderSource glShaderSource;
extern PfnGlCompileShader glCompileShader;
extern PfnGlGetShaderiv glGetShaderiv;
extern PfnGlGetShaderInfoLog glGetShaderInfoLog;
extern PfnGlDeleteShader glDeleteShader;
extern PfnGlCreateProgram glCreateProgram;
extern PfnGlAttachShader glAttachShader;
extern PfnGlLinkProgram glLinkProgram;
extern PfnGlGetProgramiv glGetProgramiv;
extern PfnGlGetProgramInfoLog glGetProgramInfoLog;
extern PfnGlUseProgram glUseProgram;
extern PfnGlDeleteProgram glDeleteProgram;
extern PfnGlGetUniformLocation glGetUniformLocation;
extern PfnGlUniformMatrix4fv glUniformMatrix4fv;
extern PfnGlUniform3f glUniform3f;
extern PfnGlUniform1f glUniform1f;

// Standard 32-bit RGBA + 24-bit depth double-buffered pixel format
// descriptor every viewer in this project uses.
PIXELFORMATDESCRIPTOR standardPixelFormatDescriptor();

// Registers a window class (if not already registered) and creates a
// window. wndProc is a parameter, not fixed to the app's own WndProc,
// specifically so callers can pass plain DefWindowProc for a throwaway
// bootstrap window -- see createGl33Context()'s comment for why that
// distinction matters.
HWND createSimpleWindow(const wchar_t* className, const wchar_t* title, int width, int height,
                         WNDPROC wndProc);

// Creates a real OpenGL 3.3 core-profile context for `hwnd` and makes it
// current, via the standard "throwaway legacy context" trick: a temporary
// window with a temporary legacy context, solely to fetch
// wglCreateContextAttribsARB (which only exists once *some* GL context is
// already current), then destroy it and create the real context on hwnd.
//
// The throwaway window deliberately uses a bare DefWindowProc, not the
// caller's real WndProc: WndProc implementations in this project call
// PostQuitMessage(0) on WM_DESTROY, and PostQuitMessage queues WM_QUIT to
// the calling *thread*, not to any specific window -- if the throwaway
// window shared the real WndProc, destroying it here would queue a
// WM_QUIT that the real window's message loop would pick up on its very
// first iteration and exit immediately, before rendering a single frame.
// (This is not a hypothetical: it is exactly the bug found and fixed
// while building the first modernized viewer, apps/viewer.)
//
// On success, also loads every GL 3.3 function pointer declared above via
// wglGetProcAddress and returns the real context (already made current on
// hOut). On failure, returns nullptr; hOut is untouched.
HGLRC createGl33Context(HWND hwnd, HDC& hOut);

// Compiles a shader from source, printing the info log to stderr on
// failure (compilation still "succeeds" in the sense that a valid-but-
// broken shader object is returned, matching plain OpenGL's own error
// model -- callers should still check linkProgram()'s own link status).
GLuint compileShader(GLenum type, const char* source);

// Links vertexShader and fragmentShader into a program, printing the info
// log to stderr on failure. Does not delete the input shaders -- callers
// own that (typically right after linking, once the program has its own
// copy).
GLuint linkProgram(GLuint vertexShader, GLuint fragmentShader);

} // namespace aether::app
