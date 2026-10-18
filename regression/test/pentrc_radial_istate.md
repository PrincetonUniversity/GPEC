# Radial LSODE failure injection

Build the fixture with the normal native compiler/dependency configuration:

```sh
cd install
make build-pentrc-radial-istate-test
```

This target builds the fixture; it does not run it. In a fresh, ordinary
PENTRC case directory containing the same input files needed by `pentrc`, run
`/path/to/gpec/bin/test_pentrc_radial_istate -2`. Exit must be nonzero and the
log must contain `PENTRC radial LSODE ISTATE=-2` and
`radial integration failed.`. It must not contain
`Radial wrapper retried failed LSODE` or
`Radial wrapper returned after a failed LSODE call`.
Repeat for -1 through -6 to cover the documented negative return statuses.

The fixture replaces only the external radial LSODE entry. Initialization,
the radial wrapper and the energy/pitch integration remain native. The
unrepaired wrapper reaches the explicit retry marker after accepting the
injected failure. Use a fresh case directory for each invocation to keep
initialization outputs separate. No scientific tolerance is changed.
