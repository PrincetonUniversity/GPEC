"""Check control-surface normals against derivatives of exported geometry.

Run after the native Solovev ideal example. Only Python's standard library
and the same NetCDF C library used by the solver are required.
"""
import argparse
import cmath
import ctypes
from ctypes.util import find_library
import math
from pathlib import Path


def read_vectors(path):
    lib = ctypes.CDLL(find_library("netcdf"))
    ncid = ctypes.c_int()

    def check(status):
        if status:
            raise RuntimeError(f"NetCDF error {status} reading {path}")

    check(lib.nc_open(str(path).encode(), 0, ctypes.byref(ncid)))
    try:
        result = []
        for name in ("R", "z", "R_n", "z_n"):
            varid, dimid, ndim = ctypes.c_int(), ctypes.c_int(), ctypes.c_int()
            size = ctypes.c_size_t()
            check(lib.nc_inq_varid(ncid, name.encode(), ctypes.byref(varid)))
            check(lib.nc_inq_varndims(ncid, varid, ctypes.byref(ndim)))
            assert ndim.value == 1, f"{name} must be a vector"
            check(lib.nc_inq_vardimid(ncid, varid, ctypes.byref(dimid)))
            check(lib.nc_inq_dimlen(ncid, dimid, ctypes.byref(size)))
            data = (ctypes.c_double * size.value)()
            check(lib.nc_get_var_double(ncid, varid, data))
            result.append(list(data))
        return result
    finally:
        check(lib.nc_close(ncid))


def derivative(values):
    """Periodic Fourier derivative, independently of bicubic derivatives."""
    n = len(values)
    # Remove the constant part before differentiating its sampled DFT.
    mean = sum(values) / n
    coefficients = []
    for k in range(1, n // 2):
        value = sum((x - mean) * cmath.exp(-2j * math.pi * k * j / n)
                    for j, x in enumerate(values)) / n
        coefficients.append((k, value))
    return [2 * sum((1j * k * value * cmath.exp(2j * math.pi * k * j / n)).real
                    for k, value in coefficients) for j in range(n)]


def check_case(path):
    r, z, rn, zn = read_vectors(path)
    assert len(r) >= 33 and len({len(x) for x in (r, z, rn, zn)}) == 1
    assert all(math.isfinite(x) for a in (r, z, rn, zn) for x in a)
    assert math.hypot(r[0] - r[-1], z[0] - z[-1]) < 1e-9
    dr, dz = derivative(r[:-1]), derivative(z[:-1])
    error = 0.0
    for a, b, nr, nz in zip(dr, dz, rn[:-1], zn[:-1]):
        length = math.hypot(a, b)
        assert length > 0
        error = max(error, math.hypot(nr - b / length, nz + a / length))
    # Differentiating the sampled fitted geometry is independent of the
    # analytic bicubic derivative and includes its discretization error.
    assert error < 1e-4, f"control-surface normal error {error:.6e}"
    print(f"PASS: control-surface normal maximum vector error {error:.6e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    check_case(args.case / "gpec_control_output_n1.nc")
