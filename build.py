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


def run(command: list[str], description: str) -> None:
    """Runs a step, echoing the command so a failure can be reproduced by
    hand rather than only through this script."""
    print(f"\n=== {description}")
    print("    " + " ".join(command), flush=True)
    result = subprocess.run(command, cwd=REPO_ROOT)
    if result.returncode != 0:
        # Exits rather than raising: a CMake or compiler failure has already
        # printed its own diagnosis, and a Python traceback on top of it
        # would bury the part that matters.
        print(f"\nFALHOU: {description} (codigo {result.returncode})", file=sys.stderr)
        sys.exit(result.returncode)


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

    run(["cmake", "-S", ".", "-B", "build"], "configurando")
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
    print("    import aether  # com python/ no PYTHONPATH, ou a partir da raiz do repo")
    return 0


if __name__ == "__main__":
    sys.exit(main())
