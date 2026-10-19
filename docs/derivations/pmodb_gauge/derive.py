"""Run the pmodb gauge identities through FortSym and record the verdicts.

Usage (FortSym checkout with a built C ABI library)::

    PYTHONPATH=$FORTSYM/python python3 docs/derivations/pmodb_gauge/derive.py

Writes ``fortsym_verdicts.json`` beside this file.  ``PROVED`` means FortSym's
native zero decision returned true; ``UNKNOWN`` and ``NONZERO`` fail an
identity and are the required outcome for a negative control.  No numerical
probe is accepted as evidence.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import fortsym  # noqa: E402
import fortsym.sympy as fs  # noqa: E402
import identities  # noqa: E402


def verdict(expr) -> str:
    value = fs.sympify(expr).is_zero
    return {True: "PROVED", False: "NONZERO"}.get(value, "UNKNOWN")


def counterexample(expr):
    """Exact rational point at which FortSym decides the residual is non-zero."""
    symbols = sorted(expr.free_symbols, key=str)
    for attempt in range(1, 6):
        point, value = {}, expr
        for j, sym in enumerate(symbols):
            num = 7 * attempt + 3 * j + 1
            point[str(sym)] = f"{num}/{5 + j}"
            value = value.subs(sym, fs.Rational(num, 5 + j))
        if value.is_zero is False:
            return point
    return None


def revision() -> str:
    root = Path(fortsym.__file__).resolve().parents[2]
    try:
        return subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main() -> int:
    records, failed = [], []
    for label, equation, residual in identities.all_identities(fs):
        start = time.perf_counter()
        v = verdict(residual)
        records.append({"label": label, "equation": equation, "verdict": v,
                        "seconds": round(time.perf_counter() - start, 3)})
        if v != "PROVED":
            failed.append(label)
        print(f"{label:6s} {v:8s} {equation}", flush=True)
    for label, equation, residual in identities.open_gaps(fs):
        v = verdict(residual)
        records.append({"label": label, "equation": equation, "verdict": v,
                        "open_gap": True})
        print(f"{label:6s} {v:8s} (open FortSym gap) {equation}", flush=True)
    for label, equation, residual in identities.negative_controls(fs):
        v = verdict(residual)
        records.append({"label": label, "equation": equation, "verdict": v,
                        "negative_control": True})
        witness = counterexample(residual) if v != "NONZERO" else None
        if witness is not None:
            records[-1]["counterexample"] = witness
        if v == "PROVED" or (v != "NONZERO" and witness is None):
            failed.append(label)
        print(f"{label:6s} {v:8s} (negative control"
              f"{', witness ' + str(witness) if witness else ''}) {equation}", flush=True)
    out = {"engine": "fortsym native zero test via C ABI (fortsym.sympy)",
           "fortsym_revision": revision(), "checks": records, "failed": failed}
    (HERE / "fortsym_verdicts.json").write_text(json.dumps(out, indent=2) + "\n")
    print("FAILED: " + " ".join(failed) if failed else "required identities proved; see explicit open-gap records")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
