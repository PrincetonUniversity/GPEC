#---------------------------------------------------------------------------
# GPEC Build System - Top-level Makefile
#
# This Makefile coordinates the build of all major GPEC modules.
# Included files:
#   DEFAULTS.inc     -- System/compiler defaults and library paths
#   SOURCES.inc      -- Source file auto-discovery
#   RULES.inc        -- Compilation rules for each module
#   TARGETS.inc      -- How to build libraries and executables
#   DEPENDENCIES.inc -- Inter-module dependencies
#
# Usage:
#   make [target]
#   Common targets:
#       all          -- Build all main executables and libraries
#       allf         -- Build all Fortran-only modules
#       clean        -- Remove built files and modules
#       clear        -- Remove example outputs (docs/examples)
#       realclean    -- Deep clean (libs, bins, modules)
#       v            -- Print build configuration
#       <module>     -- Build a specific module (e.g., pentrc, gpec)
#       print-VAR    -- Print the value of a Makefile variable (for debugging)
#
#
# See README_BUILD.md for more details.
#---------------------------------------------------------------------------


include DEFAULTS.inc
include SOURCES.inc

FFLAGS ?= -O3
CFLAGS ?= -O3

include RULES.inc

# Targets

# Top-level target
all: v neededdeps xdraw equil lsode zlange zvode orbit vacuum pentrc dcon match rdcon rmatch \
	multi sum slayer coil gpec stride

# Fortran only
allf: v neededdeps xdraw equil lsode zlange zvode orbit vacuum pentrc dcon match rdcon rmatch \
	multi sum slayer coil gpec stride

# Print verbose info
.PHONY: v
v: | FORCE
	@echo "------------------------------------------------------------"
	@echo "Defaults tailored for platform $(if $(PLATFORM),$(PLATFORM),NULL) on host $(if $(HOST),$(HOST),NULL)"
	@echo "Compiling with FC = $(FC)"
	@echo "   using math libraries $(MATHLIBS)"
	@echo "   from $(MATHDIR)"
	@echo "   including $(MATHINC)"
	@echo "   and netcdf libraries $(NETCDFLIBS) $(NETCDF_C_LIBS)"
	@echo "   from $(NETCDFDIR) $(NETCDF_C_DIR)"
	@echo "   including $(NETCDFINC)"
	@echo "   and compiler libraries $(FCLIBS)"
	@echo "   from $(FCDIR)"
	@echo "   including $(FCINC)"
	@echo "   and using OMPFLAG $(OMPFLAG)"
	@echo "   Compiling supporting modules $(NEEDED_DEPS)"
	@echo
	@echo "Set environment FC,LAPACKHOME,NETCDFHOME to change these"
	@echo "------------------------------------------------------------"

FORCE:

# Cleaning
clean:
	cd ../equil; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff
	cd ../lsode; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff
	cd ../zlange; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff
	cd ../zvode; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff
	cd ../orbit; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff orbit
	cd ../vacuum; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff
	cd ../pentrc; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff pentrc
	cd ../dcon; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff dcon
	cd ../match; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff match
	cd ../rdcon; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff rdcon
	cd ../rmatch; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff rmatch
	cd ../xdraw; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff xdraw
	cd ../sum; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff sum
	cd ../slayer; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff slayer
	cd ../coil; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff
	cd ../gpec; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff gpec
	cd ../stride; rm -f *.o *.mod *.a *.i *.original *.out *.bin sum.dat *.diff stride
	rm -rf ../bin/*


# Remove example outputs
clear:
	cd ../docs/examples/run_ideal_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/run_kinetic_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/a10_ideal_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/a10_kinetic_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/a5_tearing_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/DIIID_ideal_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/DIIID_resistive_example; rm -f *.out *.bin *.nc *.dat
	cd ../docs/examples/DIIID_kinetic_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/solovev_ideal_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/solovev_kinetic_example; rm -f *.out *.bin *.nc
	cd ../docs/examples/solovev_resistive_example; rm -f *.out *.bin *.nc *.dat

# Revert to clean state
realclean: clean clear
	rm -rf ../deps/lib/* ../deps/hdf5 ../deps/include/* ../deps/bin/* ../deps/share ../bin
	if [ -f ../deps/src/harvest/Makefile ]; then cd ../deps/src/harvest && make clean; fi
	if [ -f ../deps/src/OpenBLAS/Makefile ]; then cd ../deps/src/OpenBLAS && make clean; fi
	if [ -f ../deps/src/hdf5/Makefile ]; then cd ../deps/src/hdf5 && make clean; fi
	if [ -f ../deps/src/netcdf-c/Makefile ]; then cd ../deps/src/netcdf-c && make clean; fi
	if [ -f ../deps/src/netcdf-fortran/Makefile ]; then cd ../deps/src/netcdf-fortran && make clean; fi

# Debug
print-% : ; $(info $* is a $(flavor $*) variable set to [$($*)]) @true

# Lower level targets
include TARGETS.inc

# Dependencies
include DEPENDENCIES.inc