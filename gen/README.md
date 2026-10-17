Regenerate the supplied-rate exponential and its symbolic derivative with
FortSym. `FORTSYM_PREFIX` points to a verified FortSym build or installation
containing `include/` and `lib/`; supply its exact source revision below.

```sh
gfortran -I"$FORTSYM_PREFIX/include" gen/gen_supplied_rates.f90 \
  -L"$FORTSYM_PREFIX/lib" -lfortsym -Wl,-rpath,"$FORTSYM_PREFIX/lib" \
  -o gen_supplied_rates
./gen_supplied_rates pentrc/supplied_rates.inc FORTSYM_REVISION
```

The generated include is checked in, so normal native builds do not require
FortSym. The reader uses native spline interpolation; this generator owns only
the pointwise exponential reconstruction and its derivative. Independent
profile and invalid-input oracles run through `make check-pentrc-collisions`.
