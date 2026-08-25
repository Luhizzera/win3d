#!/usr/bin/env python3
"""One command to go from a fresh clone to a working engine.

    python build.py            configure, build, and run the test suite
    python build.py --no-test  build only
    python build.py --clean    delete build/ first, then do the above

**Why this exists rather than a line in the README.** Getting this project
running used to mean knowing three separate things: that CMake needs
configuring before building, that a multi-config generator needs
`--config Release` at *build* time rather than configure time (a
long-standing CMake trip hazard on Windows), and that the resulting
extension modules have to be found by Python afterwards. None of those are
interesting decisions a newcomer should have to make, and getting any of
them wrong produces an error that does not name its own cause. The last one
is now handled inside the package itself (see
`aether._ensure_extensions_importable`); this script handles the first two.

**What this deliberately is not**: a packaging system. It does not produce
a wheel, and `pip install aether` does not work -- that needs the C++
extensions compiled per platform and per Python version, which in turn
needs either scikit-build-core in the build chain or CI producing wheels
for each target. That is a real, separate piece of work; pretending a build
script covers it would be the kind of overstatement this project's
documentation avoids elsewhere.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
BUILD_DIR = REPO_ROOT / "build"


def run(command: list[str], description: str, capture: bool = False) -> str:
    """Runs a step, echoing the command so a failure can be reproduced by
    hand rather than only through this script.

    With `capture`, the output is both shown and returned, so a caller can
    inspect it -- used below to notice a configure that silently dropped a
    whole layer.
    """
    print(f"\n=== {description}")
    print("    " + " ".join(command), flush=True)
    if capture:
        result = subprocess.run(command, cwd=REPO_ROOT, text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        print(result.stdout, end="", flush=True)
    else:
        result = subprocess.run(command, cwd=REPO_ROOT)
    if result.returncode != 0:
        # Exits rather than raising: a CMake or compiler failure has already
        # printed its own diagnosis, and a Python traceback on top of it
        # would bury the part that matters.
        print(f"\nFALHOU: {description} (codigo {result.returncode})", file=sys.stderr)
        sys.exit(result.returncode)
    return result.stdout if capture else ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clean", action="store_true",
                        help="remove build/ before configuring")
    parser.add_argument("--no-test", action="store_true",
                        help="skip the test suite after building")
    parser.add_argument("--config", default="Release", choices=["Release", "Debug"],
                        help="build configuration (default: Release)")
    args = parser.parse_args()

    if args.clean and BUILD_DIR.exists():
        print(f"=== removendo {BUILD_DIR}")
        shutil.rmtree(BUILD_DIR)

    if shutil.which("cmake") is None:
        print("cmake nao encontrado no PATH. Instale o CMake e tente de novo.", file=sys.stderr)
        return 1

    configure_output = run(["cmake", "-S", ".", "-B", "build"], "configurando", capture=True)
    # **A configure that skips the Python bindings is not a detail to scroll
    # past.** It drops the whole orchestration layer and two of the thirteen
    # test suites -- and ctest then reports the remaining eleven passing,
    # which reads as success. `--no-tests=error` does not catch it either,
    # since tests do exist, just fewer of them. This went unnoticed once
    # already: a stale build directory carried a cached pybind11_DIR from an
    # old manual configure, so only a genuinely clean build ever hit it.
    bindings_skipped = "skipping Python bindings" in configure_output
    # --config belongs here, not on the configure line: a multi-config
    # generator picks the configuration at build time, and passing
    # CMAKE_BUILD_TYPE at configure time is silently ignored by it.
    run(["cmake", "--build", "build", "--config", args.config], f"compilando ({args.config})")

    if not args.no_test:
        # --no-tests=error rather than the default: a suite that silently
        # stops registering its tests would otherwise report success, which
        # is the one failure mode a test run must never have.
        run(["ctest", "--test-dir", "build", "-C", args.config, "--no-tests=error"],
            "rodando a suite")

    print("\n=== pronto")
    if bindings_skipped:
        print("    ATENCAO: os bindings Python foram pulados -- este build tem apenas o")
        print("             nucleo C++. O pacote 'aether' NAO vai importar, e duas suites")
        print("             de teste nao foram registradas. Instale o pybind11 no ambiente")
        print("             Python que o CMake usa, ou passe -Dpybind11_DIR=<caminho>.")
        return 1
    print("    import aether  # com python/ no PYTHONPATH, ou a partir da raiz do repo")
    return 0


if __name__ == "__main__":
    sys.exit(main())
