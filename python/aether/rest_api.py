"""Module 13 ("API": Python/C++/REST/plugins), first pass: a minimal REST
server over the `aether` package.

**Scoping decision, not asked -- unlike Module 12's dependency forks**: the
Python API (bindings for every engine layer) and the C++ API (each engine
library's own public headers) already exist and are exercised throughout
this project; REST is the one piece of Module 13 that had nothing built yet.
"Plugins" is left to Module 14, which the roadmap names as its own module
rather than folding it in here.

Built on `http.server` from the standard library -- no Flask/FastAPI/
uvicorn -- the same "don't add a dependency that isn't earning its keep"
default used everywhere else in this project (`conversational.py` made the
identical call for talking to Ollama via plain `urllib` instead of the
`requests` package). A hand-rolled request handler is a perfectly small
amount of code for the handful of routes below.

Reuses `aether.tasks` for the actual work (the same functions Module 12's
conversational layer calls) rather than reimplementing "run a cavity
simulation and turn it into JSON" a second time -- one implementation, two
front-ends, same input-validation guardrails either way.
"""

import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from aether.tasks import field_statistics, run_lid_driven_cavity

_ROUTES_POST = {
    "/simulate/lid-driven-cavity": run_lid_driven_cavity,
    "/diagnostics/field-statistics": field_statistics,
}


class _RequestHandler(BaseHTTPRequestHandler):
    def _send_json(self, status, payload):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _read_json_body(self):
        length = int(self.headers.get("Content-Length", 0))
        if length == 0:
            return {}
        return json.loads(self.rfile.read(length))

    def do_GET(self):
        if self.path == "/health":
            self._send_json(200, {"status": "ok"})
            return
        self._send_json(404, {"error": f"unknown route: GET {self.path}"})

    def do_POST(self):
        handler = _ROUTES_POST.get(self.path)
        if handler is None:
            self._send_json(404, {"error": f"unknown route: POST {self.path}"})
            return

        try:
            body = self._read_json_body()
        except json.JSONDecodeError as exc:
            self._send_json(400, {"error": f"invalid JSON body: {exc}"})
            return

        try:
            result = handler(**body)
        except (TypeError, ValueError, KeyError) as exc:
            # TypeError covers missing/unexpected fields (Python raises it
            # for a bad **body unpacking against the handler's signature),
            # ValueError covers the handler's own input validation (e.g.
            # aether.tasks.run_lid_driven_cavity's size guardrail), KeyError
            # covers a missing required key read directly out of the body.
            # All three are the caller's mistake, not a server fault, so
            # they map to 400, not 500.
            self._send_json(400, {"error": str(exc)})
            return

        self._send_json(200, result)

    def log_message(self, format, *args):
        pass  # quiet by default; BaseHTTPRequestHandler otherwise logs every request to stderr


def serve(host="127.0.0.1", port=8420):
    """Builds the server; caller decides how to run it (serve_forever() on
    the calling thread, or in a background thread for embedding/testing)."""
    return ThreadingHTTPServer((host, port), _RequestHandler)
