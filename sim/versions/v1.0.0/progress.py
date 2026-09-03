"""
progress — queryable progress for long simulations.

Every run gets its own directory `runs/<name>-<UTC timestamp>/` containing

    progress.json   status, phase, step/total, elapsed, ETA, partial results
    log.txt         append-only log
    RESULTS.md      written by the run when it finishes

`progress.json` is rewritten atomically after every tick, so
`bin/simstatus` can be called at ANY moment — while the run is going — to get
how far it is, how long is left, and the results obtained so far.

Usage
-----
    from progress import Progress

    with Progress("level1-variance", total_phases=2) as pr:
        pr.phase("closed-form", total=8)
        for k in range(8):
            ...
            pr.tick(k + 1, note=f"n={2**k}")
        pr.result("var_ratio_slope", 0.81)
        pr.write_results_md("# ...")
"""
from __future__ import annotations

import json
import os
import sys
import time
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs"


def _human(sec: float | None) -> str:
    if sec is None:
        return "?"
    sec = max(0.0, float(sec))
    if sec < 60:
        return f"{sec:.0f}s"
    if sec < 3600:
        return f"{int(sec // 60)}m{int(sec % 60):02d}s"
    return f"{int(sec // 3600)}h{int((sec % 3600) // 60):02d}m"


class Progress:
    def __init__(self, name: str, total_phases: int = 1, meta: dict | None = None):
        RUNS.mkdir(exist_ok=True)
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        self.dir = RUNS / f"{name}-{stamp}"
        self.dir.mkdir(parents=True, exist_ok=True)
        self.name = name
        self.total_phases = total_phases
        self.phase_index = 0
        self.phase_name = "(starting)"
        self.step = 0
        self.total = 0
        self.t0 = time.time()
        self.tphase = self.t0
        self.status = "running"
        self.partial: dict = {}
        self.results: dict = {}
        self.timings: list[dict] = []      # wall time of every finished phase
        self.meta = meta or {}
        # stamp the algorithm version and a code fingerprint into every run, so a
        # result can always be traced to the exact code that produced it. Best
        # effort: never let versioning break a run.
        try:
            from _algo_version import stamp as _algo_stamp
            self.meta = {**self.meta, **_algo_stamp()}
        except Exception:
            pass
        self.error: str | None = None
        self._log_fh = open(self.dir / "log.txt", "a", buffering=1)
        self.log(f"start {name} pid={os.getpid()}")
        self._flush()

    # ---------------------------------------------------------------- writing
    def _payload(self) -> dict:
        now = time.time()
        elapsed = now - self.t0
        frac = (self.step / self.total) if self.total else 0.0
        # ETA from the current phase's own rate, then scaled by phases left
        eta = None
        if self.step > 0 and self.total:
            per = (now - self.tphase) / self.step
            eta_phase = per * (self.total - self.step)
            phases_left = max(0, self.total_phases - self.phase_index)
            avg_phase = (now - self.t0) / max(1, self.phase_index)
            eta = eta_phase + phases_left * max(0.0, avg_phase)
        return {
            "name": self.name,
            "run_dir": str(self.dir),
            "status": self.status,
            "pid": os.getpid(),
            "phase_index": self.phase_index,
            "total_phases": self.total_phases,
            "phase": self.phase_name,
            "step": self.step,
            "total": self.total,
            "pct": round(100.0 * frac, 1),
            "elapsed_s": round(elapsed, 1),
            "elapsed_human": _human(elapsed),
            "eta_s": None if eta is None else round(eta, 1),
            "eta_human": _human(eta),
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(self.t0)),
            "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(now)),
            "meta": self.meta,
            "timings": self.timings,
            "phase_elapsed_s": round(now - self.tphase, 2),
            "phase_elapsed_human": _human(now - self.tphase),
            "partial": self.partial,
            "results": self.results,
            "error": self.error,
        }

    def _flush(self) -> None:
        p = self.dir / "progress.json"
        tmp = p.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(self._payload(), indent=2, default=str))
        os.replace(tmp, p)

    def log(self, msg: str) -> None:
        line = f"{time.strftime('%H:%M:%S', time.gmtime())} {msg}"
        self._log_fh.write(line + "\n")

    # ---------------------------------------------------------------- driving
    def phase(self, name: str, total: int) -> None:
        self._close_phase()
        self.phase_index += 1
        self.phase_name = name
        self.step = 0
        self.total = total
        self.tphase = time.time()
        self.log(f"phase {self.phase_index}/{self.total_phases}: {name} (total={total})")
        self._flush()

    def _close_phase(self) -> None:
        """Record the wall time of the phase that is ending."""
        if self.phase_index == 0:
            return
        el = time.time() - self.tphase
        self.timings.append({"phase": self.phase_name,
                             "index": self.phase_index,
                             "seconds": round(el, 2),
                             "human": _human(el)})
        self.log(f"phase {self.phase_index} took {_human(el)} ({el:.2f}s)")

    def tick(self, step: int | None = None, **partial) -> None:
        self.step = (self.step + 1) if step is None else step
        if partial:
            self.partial.update(partial)
        self._flush()

    def result(self, key: str, value) -> None:
        """Record a finished, quotable result. Visible immediately to simstatus."""
        self.results[key] = value
        self.log(f"result {key} = {value}")
        self._flush()

    def write_results_md(self, text: str) -> None:
        (self.dir / "RESULTS.md").write_text(text)
        self.log("wrote RESULTS.md")

    # ---------------------------------------------------------------- closing
    def timing_table_md(self) -> str:
        """Markdown timing table — every run's RESULTS.md must end with this."""
        total = time.time() - self.t0
        L = ["## Compute time", "",
             "| phase | wall time | seconds |", "|---|---|---|"]
        for t in self.timings:
            L.append(f"| {t['index']}. {t['phase']} | {t['human']} | {t['seconds']} |")
        L += [f"| **total** | **{_human(total)}** | **{total:.2f}** |", "",
              f"Machine: {os.uname().sysname} {os.uname().machine}, "
              f"python {sys.version.split()[0]}.", ""]
        return "\n".join(L)

    def done(self) -> None:
        self._close_phase()
        self.phase_index = 0            # so a later _close_phase is a no-op
        self.status = "done"
        total = time.time() - self.t0
        self.results["_total_seconds"] = round(total, 2)
        self.results["_total_human"] = _human(total)
        self.log(f"done in {_human(total)} ({total:.2f}s)")
        self._flush()

    def fail(self, msg: str) -> None:
        self.status = "failed"
        self.error = msg
        self.log(f"FAILED {msg}")
        self._flush()

    def __enter__(self) -> "Progress":
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc_type is not None:
            self.fail("".join(traceback.format_exception_only(exc_type, exc)).strip())
            self._log_fh.write("".join(traceback.format_exception(exc_type, exc, tb)))
        elif self.status == "running":
            self.done()
        self._log_fh.close()
        return False
