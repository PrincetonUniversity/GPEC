"""The LU Galerkin solve in rdcon must not depend on the BLAS thread count.

Needs a built bin/rdcon (GPEC_RDCON may point at another binary). The check
only discriminates when LAPACK/BLAS is a threaded library such as OpenBLAS;
with reference BLAS both runs agree trivially.
"""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


def delta_prime(root, exe, threads, scratch):
    """Run the a5 tearing example (solver="LU") and return the Delta_ij entry."""
    example = root / 'docs/examples/a5_tearing_example'
    with tempfile.TemporaryDirectory(dir=scratch) as directory:
        work = Path(directory)
        for name in ('equil.in', 'gfile_a5_cir', 'rdcon.in', 'vac.in'):
            shutil.copy(example / name, work)
        env = dict(os.environ, GPECHOME=str(root), OMP_NUM_THREADS=str(threads),
                   OPENBLAS_NUM_THREADS=str(threads))
        # the end of the run may fail after the matching data are written
        subprocess.run([str(exe)], cwd=work, env=env, capture_output=True)
        lines = (work / 'pest3_data_re.out').read_text().splitlines()
        return float(lines[2].split()[-1])


class RdconLuThreads(unittest.TestCase):
    def test_delta_prime_independent_of_thread_count(self):
        root = Path(__file__).resolve().parents[2]
        exe = Path(os.environ.get('GPEC_RDCON', root / 'bin/rdcon'))
        if not exe.exists():
            self.skipTest('rdcon is not built')
        scratch = Path(os.environ.get('GPEC_TEST_SCRATCH', root / 'build-test'))
        scratch.mkdir(exist_ok=True)
        values = [delta_prime(root, exe, n, scratch) for n in (1, 2, 4, 8)]
        for value in values:
            self.assertAlmostEqual(value, values[0], delta=2e-3)


if __name__ == '__main__':
    unittest.main()
