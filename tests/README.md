# Cartesian vacuum query tests

This standalone adapter builds the checkout's native libraries and the GPEC and
DCON executables, then runs the Cartesian export regression tests. It uses no
external case or private path. The upstream installation makefiles remain
unchanged.

Requirements are GNU C/Fortran, CMake 3.20 or newer, BLAS/LAPACK, NetCDF C and
Fortran development libraries, and Python 3 (standard library only). CMake's
normal prefix/library/include search paths select those dependencies. Run from
this directory with the local Fortran workflow:

```sh
FO_JOBS=2 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 fo
```

For installations without `fo`, the standalone adapter can also be run using
the standard CMake/CTest commands:

```sh
cmake -S tests -B tests/build -DCMAKE_BUILD_TYPE=Release
cmake --build tests/build --parallel 2
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ctest --test-dir tests/build --output-on-failure
```

Run the CMake commands from the checkout root. All generated inputs and native
diagnostics stay in the chosen build directory. The synthetic fixture is a
circular clockwise surface with 128 unique nodes, five poloidal modes, signed complex
normal-field coefficients and n=3. It is a bounded numerical source, not a
device equilibrium or a physical forcing convergence study.

`vacuum_query_exact` independently evaluates native CHI at double query
coordinates, then checks pickup and MSCFLD forwarding, source/stencil masks,
and identical absent/false legacy optional behavior. It includes a resolved
query inside the old snapping band and unresolved source/crossing queries.
The independent reference bypasses both coordinate construction and pickup.
The remaining tests permit finite fields and raw masked exports, and require
the specific native rejection diagnostics for each optional postprocessor.
They do not treat an arbitrary crash as a successful rejection.

These tests establish the query/output contract. They do not establish source
quadrature convergence, accuracy near a singular sheet, arbitrary malformed
derivative settings, global plasma response or physical forcing validity.
