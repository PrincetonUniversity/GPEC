#!/bin/bash -l

# This script runs all the GPEC package examples in one go.
# It should be run at the installation of a new public version so the example outputs are available and up to date.

# These ideal MHD examples are relatively fast
cd ../docs/examples/a10_ideal_example; ../../../bin/dcon; ../../../bin/stride; ../../../bin/gpec; ../../../bin/pentrc; cd ../../../install
cd ../docs/examples/solovev_ideal_example; ../../../bin/dcon; ../../../bin/stride; ../../../bin/gpec; ../../../bin/pentrc; cd ../../../install
cd ../docs/examples/DIIID_ideal_example; ../../../bin/dcon; ../../../bin/stride; ../../../bin/gpec; ../../../bin/pentrc; cd ../../../install

# These kinetic examples are slower
cd ../docs/examples/a10_kinetic_example; ../../../bin/dcon; ../../../bin/gpec; ../../../bin/pentrc; cd ../../../install
cd ../docs/examples/solovev_kinetic_example; ../../../bin/dcon; ../../../bin/gpec; ../../../bin/pentrc; cd ../../../install
cd ../docs/examples/DIIID_kinetic_example; ../../../bin/dcon; ../../../bin/gpec; ../../../bin/pentrc; cd ../../../install

# These resistive examples are a bit newer
cd ../docs/examples/DIIID_resistive_example; ../../../bin/rdcon; ../../../bin/rmatch; ../../../bin/gpec; cd ../../../install
cd ../docs/examples/solovev_resistive_example; ../../../bin/rdcon; ../../../bin/rmatch; ../../../bin/gpec; cd ../../../install
cd ../docs/examples/a5_tearing_example; ../../../bin/rdcon; ../../../bin/rmatch; cd ../../../install

# These examples have default settings preferred by developer Jong-Kyu Park
cd ../docs/examples/run_ideal_example; ../../../bin/dcon; ../../../bin/gpec; ../../../bin/pentrc; cd ../../../install
cd ../docs/examples/run_kinetic_example; ../../../bin/dcon; ../../../bin/gpec; ../../../bin/pentrc; cd ../../../install
