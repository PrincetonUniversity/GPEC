#!/usr/bin/env python3
"""
Regression check for the GPEC energy decomposition (effective_field_flag).

Run after ``dcon`` and ``gpec`` in this directory. Reads the last
GPEC_RECON block of gpec.log and checks that

  * the recon psi range covers the DCON psi range (psilow to psilim),
  * dW_p(recon-normalize) matches DCON ep(mode) within REL_TOL, and
  * the curl check reports PASS, and
  * the field and displacement tables have matching radial knots.

Exits 0 on success and 1 on failure, so it can be chained in CI.
"""

import sys

REL_TOL = 0.02  # reference run agrees to 0.4%; see README
RANGE_TOL = 1e-6


def last_recon_block(path):
    with open(path) as f:
        lines = f.read().splitlines()
    starts = [i for i, l in enumerate(lines)
              if l.startswith("GPEC_RECON final results:")]
    if not starts:
        sys.exit("FAIL: no GPEC_RECON block in " + path)
    block = []
    for l in lines[starts[-1] + 1:]:
        if l.startswith("GPEC_RECON"):
            break
        block.append(l)
    return block


def values(block, key):
    for l in block:
        if key in l:
            return [float(v) for v in l.split("=")[-1].split()]
    sys.exit("FAIL: '{}' not found in GPEC_RECON block".format(key))


def radial_modes(path):
    with open(path) as f:
        for line in f:
            if line.split()[:2] == ["psi", "m"]:
                return [(float(psi), int(mode))
                        for psi, mode, *_ in (row.split() for row in f if row.strip())]
    sys.exit("FAIL: psi and m headers not found in " + path)


def check_field_knots():
    field = radial_modes("gpec_pmodb_n1.out")
    displacement = radial_modes("gpec_xclebsch_n1.out")
    if not field or len(field) != len(displacement):
        sys.exit("FAIL: inconsistent pmodb and xclebsch table lengths")
    if any(m != n for (_, m), (_, n) in zip(field, displacement)):
        sys.exit("FAIL: inconsistent pmodb and xclebsch mode numbers")
    error = max(abs(psi - reference)
                for (psi, _), (reference, _) in zip(field, displacement))
    print("radial knot difference = {:.3e}".format(error))
    return error <= 6e-16  # Rounding of xclebsch's ES23.15 field.


def main():
    block = last_recon_block("gpec.log")
    dw = values(block, "dW_p(recon-normalize)")[0]
    ep = values(block, "DCON ep(")[0]
    rlo, rhi = values(block, "recon psi range")
    dlo, dhi = values(block, "DCON psi range")
    rel = abs(dw - ep) / abs(ep)
    verdict = [l.split(":")[-1].strip() for l in block if "C verify" in l]
    range_ok = abs(rlo - dlo) <= RANGE_TOL and abs(rhi - dhi) <= RANGE_TOL
    print("recon psi range       = {:.6e} {:.6e}".format(rlo, rhi))
    print("DCON psi range        = {:.6e} {:.6e}".format(dlo, dhi))
    print("dW_p(recon-normalize) = {: .8e}".format(dw))
    print("DCON ep               = {: .8e}".format(ep))
    print("relative difference   = {:.3e} (tol {:.1e})".format(rel, REL_TOL))
    print("curl check            = {}".format(verdict[-1] if verdict
                                              else "not run"))
    ok = (range_ok and rel <= REL_TOL and bool(verdict)
          and verdict[-1] == "PASS")
    knots_ok = check_field_knots()
    print("PASS" if ok and knots_ok else "FAIL")
    return 0 if ok and knots_ok else 1


if __name__ == "__main__":
    sys.exit(main())
