"""Module 12's "conversational" sub-area (first pass): a small tool-using
chat layer over the `aether` package, talking to a *local* Ollama server
rather than a cloud LLM API.

Chosen deliberately over Anthropic/OpenAI-style cloud APIs after the user
was asked explicitly (this is a genuine fork, unlike most dependency
choices in this project): a local model on the project's own RTX 5080
(already set up for CUDA in Module 10) needs no API key, costs nothing per
call, and sends no data off the machine -- the closest available option to
this project's standing "no external dependency" default, applied to a
sub-area that structurally cannot be fully dependency-free (some LLM has to
generate the natural-language side of "conversational").

Talks to Ollama's REST API (http://localhost:11434 by default) using only
`urllib` from the standard library -- no `ollama`/`openai` pip package --
the same "don't add a dependency that isn't earning its keep" instinct
used everywhere else in this project.

**The LLM never computes a physical quantity itself.** It only ever decides
*which* aether function to call and with *what* arguments; the actual
number always comes from actually running the real solver/diagnostic in
this process, not from the model's own text. `AetherAssistant.ask()`
returns both the model's natural-language answer and the raw tool-call
results (tool name, arguments, and the real computed result) precisely so
a caller can check the reported number against the engine's own output
instead of trusting the model's prose -- the same "don't just trust
plausible-looking output" discipline this project applies to its own
numerical code.

Deliberately NOT imported by `aether/__init__.py`: unlike the GPU module
(whose availability can be checked at import time via whether the compiled
extension exists), there is no way to know whether a local Ollama server is
actually running without trying a network call, so this stays an explicit,
opt-in import (`from aether.conversational import AetherAssistant`) rather
than something that could make importing the base `aether` package itself
depend on an external process being up.
"""

import json
import urllib.error
import urllib.request

from aether.tasks import run_lid_driven_cavity


class OllamaClient:
    """Thin wrapper around Ollama's /api/chat endpoint."""

    def __init__(self, model="llama3.2", host="http://localhost:11434"):
        self.model = model
        self.host = host

    def chat(self, messages, tools=None):
        payload = {"model": self.model, "messages": messages, "stream": False}
        if tools:
            payload["tools"] = tools
        data = json.dumps(payload).encode("utf-8")
        request = urllib.request.Request(
            f"{self.host}/api/chat", data=data, headers={"Content-Type": "application/json"}
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                body = json.loads(response.read().decode("utf-8"))
        except urllib.error.URLError as exc:
            raise RuntimeError(
                f"Could not reach the local Ollama server at {self.host} -- is 'ollama serve' running?"
            ) from exc
        return body["message"]


# Each entry pairs a JSON-schema tool description (what the model sees) with
# the real Python callable that actually computes the answer (what runs).
# The callable itself -- including its size guardrail, see
# aether.tasks.MAX_CAVITY_CELL_STEPS's own comment for why that exists --
# lives in aether.tasks, shared with Module 13's REST API rather than
# duplicated here.
_TOOLS = {
    "run_lid_driven_cavity": {
        "callable": run_lid_driven_cavity,
        "schema": {
            "type": "function",
            "function": {
                "name": "run_lid_driven_cavity",
                "description": (
                    "Runs a real 2D lid-driven cavity incompressible Navier-Stokes "
                    "simulation (Aether's own solver -- not an estimate or a lookup) "
                    "and reports the maximum mass-conservation divergence, the "
                    "simulated time reached, and the pressure field's checkerboard "
                    "(odd-even decoupling) index after the given number of explicit "
                    "time steps."
                ),
                "parameters": {
                    "type": "object",
                    "properties": {
                        "nx": {"type": "integer", "description": "grid cells along x"},
                        "ny": {"type": "integer", "description": "grid cells along y"},
                        "reynolds_number": {"type": "number", "description": "lid_velocity * length / viscosity"},
                        "steps": {"type": "integer", "description": "number of explicit time steps to run"},
                    },
                    "required": ["nx", "ny", "reynolds_number", "steps"],
                },
            },
        },
    },
}


class AetherAssistant:
    def __init__(self, model="llama3.2", host="http://localhost:11434"):
        self.client = OllamaClient(model=model, host=host)

    def ask(self, question):
        tool_schemas = [tool["schema"] for tool in _TOOLS.values()]
        messages = [{"role": "user", "content": question}]

        message = self.client.chat(messages, tools=tool_schemas)
        tool_calls = message.get("tool_calls") or []

        executed = []
        if tool_calls:
            messages.append(message)
            for call in tool_calls:
                name = call["function"]["name"]
                arguments = call["function"]["arguments"]
                if name not in _TOOLS:
                    raise RuntimeError(f"model requested unknown tool '{name}'")
                try:
                    result = _TOOLS[name]["callable"](**arguments)
                except (ValueError, TypeError) as exc:
                    # A guardrail (or a malformed argument) rejected the
                    # call -- reported back to the model as a tool result,
                    # same as a real answer, rather than raised out of
                    # ask() and left to crash or hang the caller. See
                    # aether.tasks.MAX_CAVITY_CELL_STEPS's comment for why
                    # this matters in practice, not just in theory.
                    result = {"error": str(exc)}
                executed.append({"name": name, "arguments": dict(arguments), "result": result})
                messages.append({"role": "tool", "content": json.dumps(result)})

            message = self.client.chat(messages, tools=tool_schemas)

        return {"answer": message.get("content", ""), "tool_calls": executed}
