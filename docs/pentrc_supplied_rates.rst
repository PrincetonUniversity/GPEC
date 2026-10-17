PENTRC supplied collision rates
==============================

``kinetic_format = "rates-v1"`` in ``PENT_INPUT`` selects a nine-column
``kinetic_file``. The default ``"legacy"`` retains the six-column format and
its existing main-ion/impurity collision approximation. The columns are::

    psi_n  n_i  n_e  T_i  T_e  OmegaE  logLambda  nu_i  nu_e

Densities are in m^-3, temperatures in eV, physical electric rotation
``OmegaE`` in rad/s, the Coulomb logarithm dimensionless, and the supplied
base Krook damping rates in s^-1. No 2*pi conversion is applied. All densities,
temperatures, logarithms and rates must be positive and finite. A selected-ion
charge density greater than the total electron density is rejected. Exact
zero rotation is retained in this format. The coordinate orientation is
unchanged: provide ``OmegaE`` in native axes with ``wefac=1``, or declare the
adapter's common-to-native axis conversion through ``wefac``. Supplied-rate
mode requires ``wpfac=1``.

The supplied rates bypass the legacy impurity-rate formula at both native
query interfaces. Native splines interpolate their logarithms on the input
radial grid, and the generated exponential kernel restores their values and
derivatives. This preserves positivity between positive input samples;
overflow, underflow and nonfinite results are rejected. Primitive drive
profiles retain the existing interpolation. As in legacy mode, ``nfac`` and
``tfac`` alter drive profiles only, leaving supplied collision rates unchanged.

One native call calculates one population. A common D/T/e mixture therefore
uses an adapter that emits three calls. For D and T, column ``n_i,T_i`` holds
that isotope's density and temperature, and ``mi`` is respectively 2 or 3
in proton-mass units. Set ``zi=zimp=1``; ``mimp`` may identify the other isotope.
For electrons set ``electron=.true.`` and use the same total ``n_e,T_e``.
All three calls carry the same physical ``OmegaE`` and common ``logLambda``.
Supply each population's rate from the declared full background mixture,
including its densities, temperatures, charges and masses. Keep that mixture
and rate convention in the adapter's provenance record: this profile contains
the selected drive and scalar rates, not the full background distribution.
This is a diagonal scalar response, not a coupled Fokker--Planck solve.

The existing pitch and energy models are unchanged. Passing particles use
the supplied base rate; trapped particles use the base rate divided by
``2*epsilon``. ``nutype="krook"`` then uses this rate without energy
dependence. ``nutype="harmonic"`` additionally multiplies it by
``(1+ell_eff**2/4)*x**(-3/2)``, with ``x=E/T``, integer ``ell_eff=ell`` for
trapped particles and ``ell_eff=ell+n*q`` for passing particles. ``nufac``
multiplies the resulting width. Do not include these geometry, harmonic,
energy or ``nufac`` factors in the supplied base columns. An energy-only
``x**(-3/2)`` operator and a geometry-uniform constant operator are not
provided by this interface. ``nutype="zero"`` selects the collisionless model;
``"small"`` retains its numerical width proportional to ``abs(OmegaE)``.

``pentrc_kinetics.out`` records the format when ``indebug=.true.``. The native
NetCDF output records ``kinetic_format``, the consumed ``collision_model`` and
``collision_multiplier``, along with the existing selected mass, charge,
densities, temperatures, logarithm, base rates and trapped ``nu_eff`` profile.

Run the focused native reader/query regressions with::

    cd install
    make check-pentrc-collisions

The supplied-rate oracle checks unequal D/T drive profiles, both ion and
electron rates and their derivatives, physical zero rotation, and invalid
input rejection in both Debug and ``-O3 -ffast-math`` builds. The generated
rate reconstruction is reproducible from ``gen/gen_supplied_rates.f90`` using
FortSym; the checked-in ``pentrc/supplied_rates.inc`` is used by native builds.
