"""Replay a local rational-surface window and check its overlaps directly.

Run after the native Solovev ideal DCON example. The oracle contracts the
exported, area-weighted right singular vectors with the external field;
it does not use GPEC's in-memory overlap vectors.
"""
import argparse
import cmath
from pathlib import Path
import re
import shutil
import subprocess


def check_case(case):
    local = (case / "gpec_singcoup_svd_local_n1.out").read_text()
    global_svd = (case / "gpec_singcoup_svd_n1.out").read_text()
    count = lambda text: int(re.search(r"msing\s*=\s*(\d+)", text)[1])
    nlocal = count(local)
    assert 0 < nlocal < count(global_svd), "A strict local subset is required"

    control = (case / "gpec_control_n1.out").read_text()
    # The output-coordinate table is appended only when a conversion is
    # needed. Otherwise the first, working-coordinate table is the output.
    if "jac_out =" in control:
        control = control.split("jac_out =", 1)[1]
    control = control.split("real(bin)", 1)[1].splitlines()[1:]
    field = {}
    for line in control:
        words = line.split()
        if not words:
            if field:
                break
            continue
        field[int(words[0])] = complex(float(words[1]), float(words[2]))

    vectors = []
    for model in local.split("Right singular vectors")[1:]:
        modes = []
        for mode in model.split("mode =")[1:]:
            vector = {}
            for line in mode.splitlines()[1:]:
                words = line.split()
                if len(words) == 3:
                    try:
                        vector[int(words[0])] = complex(
                            float(words[1]), float(words[2]))
                    except ValueError:
                        pass
            assert vector, "Missing local singular vector"
            modes.append(vector)
        assert len(modes) == nlocal
        vectors.append(modes)
    assert len(vectors) == 5, "All five coupling models must be exercised"

    result = (case / "gpec_singfld_n1.out").read_text()
    result = result.split("Local overlap fields,", 1)[1]
    rows = []
    for line in result.splitlines():
        words = line.split()
        if len(words) == 16 and words[0].isdigit():
            rows.append([float(word) for word in words[1:]])
    assert len(rows) == nlocal, "Missing local overlap modes"
    max_error = 0.0
    for model, modes in enumerate(vectors):
        for index, vector in enumerate(modes):
            expected = sum(value.conjugate() * field[m]
                           for m, value in vector.items())
            actual = complex(*rows[index][3 * model:3 * model + 2])
            assert cmath.isfinite(actual) and cmath.isfinite(expected)
            error = abs(actual - expected) / max(abs(expected), 1.0e-12)
            max_error = max(max_error, error)
            # ASCII vectors and fields have nine significant digits.
            assert error < 2.0e-6, (model, index, actual, expected, error)
    print(f"PASS: {nlocal} local modes, five models; "
          f"maximum direct overlap relative error {max_error:.3e}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-case", type=Path)
    parser.add_argument("--work-dir", type=Path)
    parser.add_argument("--binary", type=Path)
    args = parser.parse_args()
    if args.check_case:
        check_case(args.check_case.resolve())
        return
    if not args.work_dir:
        parser.error("--work-dir is required for the Solovev replay")
    root = Path(__file__).resolve().parents[2]
    binary = args.binary.resolve() if args.binary else root / "bin/gpec"
    fixture = root / "docs/examples/regression_solovev_ideal_example"
    for name in ("euler.bin", "psi_in.bin"):
        assert (fixture / name).is_file(), "Run the native Solovev DCON first"
    for enabled in (True, False):
        case = args.work_dir.resolve() / ("singfld_on" if enabled else "singfld_off")
        case.mkdir(parents=True)
        for name in ("euler.bin", "psi_in.bin"):
            shutil.copy2(fixture / name, case / name)
        for pattern in ("*.in", "*.txt", "*.dat"):
            for source in fixture.glob(pattern):
                shutil.copy2(source, case / source.name)
        deck = (case / "gpec.in").read_text()
        deck = re.sub(r"dcon_dir\s*=\s*'[^']*'",
                      "dcon_dir = ''", deck, flags=re.I)
        deck = re.sub(r"&GPEC_OUTPUT", "&GPEC_OUTPUT\n   opsi2=0.5",
                      deck, flags=re.I)
        if not enabled:
            deck = re.sub(r"\bsingfld_flag\s*=\s*t\b", "singfld_flag=f",
                          deck, flags=re.I)
        (case / "gpec.in").write_text(deck)
        with (case / "gpec.log").open("w") as log:
            subprocess.run([binary], cwd=case, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        assert "GPEC STOP => Normal termination." in (case / "gpec.log").read_text()
        if enabled:
            check_case(case)
        else:
            assert not (case / "gpec_singfld_n1.out").exists()
            print("PASS: local coupling with singular-field output disabled")


if __name__ == "__main__":
    main()
