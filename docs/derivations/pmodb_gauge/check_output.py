"""Numeric oracle for the pmodb gauge identities on the TC24 n=3 GPEC run.

Usage: PYTHONPATH=$RMP_TORQUE/src python3 check_output.py RUN_DIR

Requires numpy, scipy, and the gpec_ascii readers from the referenced
rmp_torque Boozer-gauge-total branch. Numerical evidence only; see derive.py
for the native FortSym proofs.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from rmp_torque.gpec_ascii import (read_pmodb, read_xclebsch, read_equilibrium_grid,
                                   reconstruct_advective_spectrum)

R = Path(sys.argv[1])
MU0 = 4e-7 * np.pi


def table(path, ncol):
    rows = []
    for line in open(path):
        f = line.split()
        if len(f) != ncol:
            continue
        try:
            rows.append([float(x) for x in f])
        except ValueError:
            pass
    return np.array(rows)


d = table(R / "gpec_diagnostics_pmodb_n3.out", 16)
psi = np.unique(d[:, 0]); ms = np.unique(d[:, 1]).astype(int)
shape = (psi.size, ms.size)
c = lambda k: (d[:, k] + 1j * d[:, k + 1]).reshape(shape)
lagb, llagb, lll, bdiv, bk = c(2), c(4), c(6), c(8), c(12)
sp = read_pmodb(R / "gpec_pmodb_n3.out")
xc = read_xclebsch(R / "gpec_xclebsch_n3.out")
eq = read_equilibrium_grid(R / "idcon_equil.out")
q = np.interp(psi, eq.psi_poloidal, eq.q)
nres = 3 * q
far = np.abs(nres - np.round(nres)) > 0.15
rel = lambda a, b: np.linalg.norm(a - b, axis=1) / np.linalg.norm(b, axis=1)
bands = [(0.0, 0.8), (0.8, 0.9), (0.9, 0.98)]


def report(name, a, b):
    r = rel(a, b)
    parts = []
    for lo, hi in bands:
        sel = (psi >= lo) & (psi < hi)
        parts.append(f"{lo:.1f}-{hi:.2f}: med {np.median(r[sel]):.4f} "
                     f"(far {np.median(r[sel & far]):.4f})")
    print(f"{name:34s}", " | ".join(parts))


print("diag grid", shape, "pmodb grid", sp.lagrangian.shape)
assert np.allclose(sp.psi_poloidal, psi)
report("pmodb lagb vs -(Bdiv+Bk) [cols]", sp.lagrangian,
       -(sp.b_div_xi_perp + sp.b_kappa_xi_perp))
report("diag lagb vs pmodb lagb", lagb, sp.lagrangian)
report("diag lagb vs lllagb=-(Bdiv+Bk)", lagb, lll)
report("diag llagb vs lllagb", llagb, lll)
report("gauge term (lagb-llagb) / llagb", lagb - llagb + llagb, llagb)

# Finding (ii): the pressure-coefficient fit.
fl = []
sec = False
for line in open(R / "idcon_equil.out"):
    if "Flux functions:" in line:
        sec = True
        continue
    f = line.split()
    if sec and len(f) == 6:
        try:
            fl.append([float(x) for x in f])
        except ValueError:
            pass
fl = np.array(fl)
dp = CubicSpline(fl[:, 0], fl[:, 2]).derivative()(xc.psi_poloidal)
th = eq.theta_radians
ph = np.exp(1j * np.outer(xc.modes, th))
bs = CubicSpline(eq.psi_poloidal, eq.b_magnitude, axis=0)(xc.psi_poloidal)


def spec(x):
    r = np.fft.fft((x @ ph) / bs, axis=1) / len(th)
    return np.stack([r[:, int(m) % len(th)] for m in xc.modes], axis=1)


adv = reconstruct_advective_spectrum(xc, eq)        # xi_g . grad B (xi^zeta = 0)
A1 = spec(adv)
A2_wrong = spec(MU0 * xc.xi_psi * dp[:, None])      # rmp_torque #56 fit: missing 1/B
_r = np.fft.fft(MU0 * dp[:, None] * (xc.xi_psi @ ph) / bs**2, axis=1) / len(th)
A2 = np.stack([_r[:, int(m) % len(th)] for m in xc.modes], axis=1)   # (mu0 p' xi^psi/B)/B
A3 = spec(llagb - lagb)                             # (xi_perp^zeta/q) d_theta B
y = spec(sp.b_kappa_xi_perp)
sel = far[:, None] & (np.abs(xc.modes) <= 12)[None, :]
for name, cols in (("adv, p/B [#56]", (A1, A2_wrong)), ("adv, p", (A1, A2)), ("adv, p, gauge", (A1, A2, A3))):
    M = np.stack([a[sel] for a in cols], axis=1)
    coef, *_ = np.linalg.lstsq(M, y[sel], rcond=None)
    res = np.linalg.norm(M @ coef - y[sel]) / np.linalg.norm(y[sel])
    print(f"fit Bkxprp ~ {name:14s}: coef", np.round(coef, 4), f"resid {res:.4f}")
M = A1 + A2 + A3
print("fixed coefficients (1,1,1): resid %.4f" % (
    np.linalg.norm(M[sel] - y[sel]) / np.linalg.norm(y[sel])))

print("=== consistency of the regressors (far surfaces, |m|<=12)")
eul = sp.eulerian
def r2(a, b):
    return np.linalg.norm(a[sel] - b[sel]) / np.linalg.norm(b[sel])
print("adv (xclebsch) vs lagb-eulb      %.4f" % r2(A1, spec(sp.lagrangian - eul)))
print("p-term vs Bk-(llagb-eulb)        %.4f" % r2(A2, spec(sp.b_kappa_xi_perp - llagb + eul)))
print("Bk vs (llagb-eulb)+p-term        %.4f" % r2(y, spec(llagb - eul) + A2))

T = spec(sp.b_kappa_xi_perp - llagb + eul)
scale = np.vdot(A2[sel], T[sel]) / np.vdot(A2[sel], A2[sel])
print("pressure coefficient against corrected regressor:", scale)
