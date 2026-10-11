"""Endpoint extrapolation must preserve exact low-degree profiles."""
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


class SplineEndpointPolynomials(unittest.TestCase):
    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is required')
    def test_clustered_grid_preserves_polynomial_derivatives_and_flux(self):
        root = Path(__file__).resolve().parents[2]
        scratch = Path(os.environ.get('GPEC_TEST_SCRATCH', root / 'build-test'))
        scratch.mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=scratch) as directory:
            work = Path(directory)
            for name in ('io', 'local', 'utils', 'spline', 'cspline'):
                source = root / 'equil' / (name + '.f')
                if name == 'spline' and os.environ.get('GPEC_SPLINE_SOURCE'):
                    source = Path(os.environ['GPEC_SPLINE_SOURCE'])
                subprocess.run(['gfortran', '-O2', '-ffixed-line-length-none',
                                '-c', str(source)], cwd=work, check=True,
                               capture_output=True)
            driver = root / 'regression/test/spline_endpoint_polynomials.f90'
            subprocess.run(['gfortran', '-O2', str(driver),
                            *map(str, sorted(work.glob('*.o'))),
                            '-llapack', '-lblas', '-o', str(work / 'probe')],
                           cwd=work, check=True, capture_output=True)
            result = subprocess.run([str(work / 'probe')], cwd=work, check=True,
                                    capture_output=True, text=True)
            for line in result.stdout.splitlines():
                if line.startswith('complex'):
                    _, resolution, derivative = line.split()
                    with self.subTest(resolution=resolution, profile='complex'):
                        self.assertLess(float(derivative), 1e-7)
                    continue
                resolution, degree, derivative, endpoint, integral = line.split()
                with self.subTest(resolution=resolution, degree=degree):
                    self.assertLess(float(derivative), 1e-7)
                    self.assertLess(float(endpoint), 1e-8)
                    self.assertLess(float(integral), 1e-10)


if __name__ == '__main__':
    unittest.main()
