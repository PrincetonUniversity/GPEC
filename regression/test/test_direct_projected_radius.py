"""Check that the delivered direct chart stays on its input flux surface."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class DirectProjectedRadius(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is required')
    def test_chart_radius_follows_input_flux_surface(self):
        root = Path(__file__).resolve().parents[2]
        scratch = Path(os.environ.get('GPEC_TEST_SCRATCH', root / 'build-test'))
        scratch.mkdir(exist_ok=True)
        names = ('io', 'local', 'utils', 'spline', 'bicube', 'cspline',
                 'global', 'grid', 'inverse', 'direct', 'read_eq')
        sources = [root / 'equil' / (name + ('.f90' if name == 'grid' else '.f'))
                   for name in names]
        direct_override = os.environ.get('GPEC_DIRECT_SOURCE')
        if direct_override:
            sources[names.index('direct')] = Path(direct_override)
        sources += sorted((root / 'lsode').glob('*.f'))
        fixture = root / 'docs/examples/a5_tearing_example/gfile_a5_cir'
        driver = root / 'regression/test/direct_projected_radius.f90'
        with tempfile.TemporaryDirectory(dir=scratch) as directory:
            work = Path(directory)
            for source in sources:
                subprocess.run(['gfortran', '-O2', '-fallow-argument-mismatch',
                                '-frecursive', '-ffixed-line-length-none', '-c',
                                str(source)], cwd=work, check=True,
                               capture_output=True)
            subprocess.run(['gfortran', '-O2', str(driver),
                            *map(str, sorted(work.glob('*.o'))),
                            '-llapack', '-lblas', '-o', str(work / 'probe')],
                           cwd=work, check=True, capture_output=True)
            result = subprocess.run([str(work / 'probe'), str(fixture)],
                                    cwd=work, check=True, capture_output=True,
                                    text=True)
            residual = float(result.stdout.splitlines()[-1])
            self.assertLess(residual, 1e-5)


if __name__ == '__main__':
    unittest.main()
