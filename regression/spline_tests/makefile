# Unit tests for equil/spline.f. Build GPEC first (install/make), then: make run
FC ?= gfortran
GPEC = ../..
LIBS ?= -L$(LAPACKHOME)/lib -llapack -lblas

test_spline_pchip: test_spline_pchip.f90 $(GPEC)/deps/lib/libequil.a
	$(FC) $(FFLAGS) -I$(GPEC)/equil -o $@ $< $(GPEC)/deps/lib/libequil.a $(LIBS)

run: test_spline_pchip
	./test_spline_pchip

clean:
	rm -f test_spline_pchip *.o *.mod
