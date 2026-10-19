"""Cross-check the native gauge test's Hamada->Boozer result with rmp_torque.

Run pmodb_target_gauge first. Set PYTHONPATH to rmp_torque #55 source.
"""
import numpy as np
from rmp_torque.boozer_gauge import Spectrum, hamada_to_boozer

rows = np.loadtxt('pmodb_hamada_to_boozer.txt')
modes = np.arange(-12, 13)
values = np.zeros(25, complex)
values[13] = 0.2
bvalues = np.zeros(25, complex)
bvalues[12] = 3
bvalues[11] = bvalues[13] = -0.3
# GPEC exp(i(m theta-n zeta)) maps to ITP n=-2; both angles here in radians.
expected = hamada_to_boozer(Spectrum(modes, values), Spectrum(modes,bvalues), -2, 1/1.3, 12)
actual = rows[:,2]+1j*rows[:,3]
reference = expected.evaluate(2*np.pi*rows[:,0])
np.testing.assert_allclose(actual, reference, atol=2e-11, rtol=2e-11)
print('PASS: native Hamada->Boozer output agrees with rmp_torque #55')
