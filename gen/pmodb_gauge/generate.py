"""Prove and emit the inverse field-line coefficient (angles in turns)."""
from pathlib import Path
import fortsym
import fortsym.sympy as fs
s = fs.Symbol("s")
weight = 1/(2*fs.pi*s)
assert (2*fs.pi*s*weight-1).is_zero is True
source = fortsym.emit_kernel("pmodb_gauge_inverse", ["s"], {"weight":weight},
                             generator="gen/pmodb_gauge/generate.py")
root = Path(__file__).resolve().parents[2]
(root/"gpec/pmodb_gauge_kernel.f90").write_text(
    "module pmodb_gauge_kernel\ncontains\n"+source+"\nend module pmodb_gauge_kernel\n")
