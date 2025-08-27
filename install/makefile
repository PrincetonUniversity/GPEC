include DEFAULTS.inc
include SOURCES.inc

FFLAGS ?= -O3
CFLAGS ?= -O3

include RULES.inc

# Targets

# Top-level target
all: v harvest equil lsode zlange zvode orbit vacuum pentrc dcon match rdcon rmatch \
	multi sum xdraw slayer coil gpec stride

# Fortran only
allf: v harvest equil lsode zlange zvode orbit vacuum pentrc dcon match rdcon rmatch \
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
	@echo
	@echo "Set environment FC,LAPACKHOME,NETCDFHOME to change these"
	@echo "------------------------------------------------------------"

FORCE:

# Cleaning
clean:
	cd ../equil; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff
	cd ../lsode; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff
	cd ../zlange; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff
	cd ../zvode; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff
	cd ../orbit; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff orbit
	cd ../vacuum; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff
	cd ../dcon; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff dcon
	cd ../match; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff match
	cd ../rdcon; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff rdcon
	cd ../rmatch; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff rmatch
	cd ../xdraw; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff xdraw
	cd ../sum; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff sum
	cd ../slayer; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff slayer
	cd ../coil; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff
	cd ../gpec; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff gpec
	cd ../stride; rm -f *.o *.mod *.a *.original *.out *.bin sum.dat *.diff stride
	rm -rf ../bin/*
	cd ../harvest; make clean

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
	rm -rf ../lib ../bin

# Debug
print-% : ; $(info $* is a $(flavor $*) variable set to [$($*)]) @true

# Lower level targets
include TARGETS.inc

# Dependencies
include DEPENDENCIES.inc