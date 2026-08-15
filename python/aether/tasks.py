"""Operations exposed to *external* interfaces -- Module 12's conversational
layer and Module 13's REST API both wrap the same underlying calls here,
instead of each keeping its own copy of "how to run a cavity simulation and
turn it into a JSON-able result". One implementation, one guardrail, two
front-ends.

Each function here does the same thing a hand-written script already could
via the `aether` package directly: build a real solver, run it, read off
real diagnostics. Nothing is estimated or looked up -- these exist purely to
give an external caller (an LLM choosing tool arguments, an HTTP client
posting JSON) a single, input-validated entry point.
"""

import aether

# Guardrail found necessary by actually testing the conversational layer
# against a local LLM, not anticipated in advance: asked a purely
# conceptual question ("what is a Reynolds number?"), llama3.2 (a small
# 3B model) still decided to call this operation anyway, picking
# nx=128, ny=64, steps=10000 out of nowhere -- an 8192-cell grid run for
# 10000 explicit steps, a real, legitimately expensive simulation
# (confirmed: it pegged a CPU core for minutes before being killed). A
# caller -- whether an LLM guessing at arguments or an HTTP client with a
# typo -- cannot be trusted to only request reasonably-sized work, so this
# module itself, not whoever calls it, refuses anything above a sane
# interactive bound before running it.
MAX_CAVITY_CELL_STEPS = 2_000_000  # nx*ny*steps; ~10x the 20x20x500 example this module was validated against


def run_lid_driven_cavity(nx, ny, reynolds_number, steps, lid_velocity=1.0, length=1.0):
    """Runs the real LidDrivenCavitySolver2D and reports diagnostics computed
    by the engine itself (max_divergence from the solver, checkerboard_index
    from Module 12.2's FlowDiagnostics). Raises ValueError on invalid or
    oversized input -- callers are expected to turn that into whatever error
    shape fits their own protocol (a REST 400, a tool-call error message)."""
    nx, ny, steps = int(nx), int(ny), int(steps)
    reynolds_number = float(reynolds_number)
    lid_velocity = float(lid_velocity)
    length = float(length)

    if nx <= 0 or ny <= 0 or steps <= 0:
        raise ValueError("nx, ny and steps must all be positive")
    if nx * ny * steps > MAX_CAVITY_CELL_STEPS:
        raise ValueError(
            f"requested simulation is too large (nx*ny*steps={nx * ny * steps} exceeds the "
            f"{MAX_CAVITY_CELL_STEPS} guardrail) -- ask for a smaller grid or fewer steps"
        )

    viscosity = lid_velocity * length / reynolds_number

    solver = aether.LidDrivenCavitySolver2D(nx, ny, length, length, viscosity, lid_velocity)
    for _ in range(steps):
        solver.step(solver.stable_time_step())

    pressure = [solver.pressure(i, j) for j in range(ny) for i in range(nx)]
    return {
        "max_divergence": solver.max_divergence(),
        "simulated_time": solver.time(),
        "checkerboard_index": aether.checkerboard_index(pressure, nx, ny),
    }


def field_statistics(field):
    """Wraps FlowDiagnostics' computeStatistics for external callers that
    only have a plain list of numbers, not a live solver object."""
    stats = aether.compute_statistics(field)
    return {"min_value": stats.min_value, "max_value": stats.max_value, "mean": stats.mean}
