"""Algorithm versioning for the simulation code.

Every run records which version of the algorithms produced it, and a fingerprint
of the exact code, so that a result in `runs/` can always be traced back to the
code that made it.  `progress.Progress` stamps both into every run's meta
automatically, so no runner has to remember to do it.

Bump `ALGO_VERSION` (and snapshot the code) with `bin/newalgo`; the manifest is
`ALGORITHMS.md`.  Semantics mirror the paper's semver:

    patch  bug fix or refactor with no change to any numerical result
    minor  a new algorithm, or a change that alters results in an intended way
    major  a reversed result or a construction replaced

The fingerprint is a short SHA over all `sim/*.py` and `sim/cpp/*.cpp` (their
paths and contents), so any edit to any algorithm changes it.  It is the machine
check; `ALGO_VERSION` is the human label.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

ALGO_VERSION = "1.3.1"

_SIM = Path(__file__).resolve().parent


def code_fingerprint() -> str:
    """Short SHA over every algorithm source file (path + bytes)."""
    h = hashlib.sha256()
    files = sorted(_SIM.glob("*.py")) + sorted((_SIM / "cpp").glob("*.cpp"))
    for f in files:
        try:
            h.update(f.relative_to(_SIM).as_posix().encode())
            h.update(f.read_bytes())
        except OSError:
            continue
    return h.hexdigest()[:12]


def stamp() -> dict:
    """The dict merged into every run's meta by progress.Progress."""
    return {"algo_version": ALGO_VERSION, "code_sha": code_fingerprint()}
