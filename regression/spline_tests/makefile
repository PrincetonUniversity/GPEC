# Unit tests for equil/spline.f. Build GPEC first (install/make), then: make run
FC ?= gfortran
GPEC = ../..
LIBS ?= -L$(LAPACKHOME)/lib -llapack -lblas

test_spline_hermite: test_spline_hermite.f90 $(GPEC)/deps/lib/libequil.a
	$(FC) $(FFLAGS) -I$(GPEC)/equil -o $@ $< $(GPEC)/deps/lib/libequil.a $(LIBS)

run: test_spline_hermite
	./test_spline_hermite

clean:
	rm -f test_spline_hermite *.o *.mod
