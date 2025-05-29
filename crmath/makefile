# Makefile for Fortran correctly rounded math library

include ../install/DEFAULTS.inc
# If using ifort, the proper syntax is -recursive.
ifneq (,$(findstring ifort,$(FC)))
	STANDARDFLAG=-stand=f08
else ifneq (,$(findstring ifx,$(FC)))
	STANDARDFLAG=-stand=f08
else
	STANDARDFLAG=-std=f2008
endif

IFLAGS = -I$(CRLIBMINC)

OBJS = \
	crmath.o \

F90 = $(FC) -O2 -fPIC $(STANDARDFLAG) $(IFLAGS)

%.o: %.f90
	$(F90) -c $<

%.o: %.F90
	$(F90) -c $<

all : libcrmath

clean :
	rm -f $(OBJS) *.mod *.smod *.a

libcrmath : $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libcrmath.a $(OBJS)

# dependencies

crmath_r_sp.o : crmath.mod
crmath_c_sp.o : crmath.mod
crmath_r_dp.o : crmath.mod
crmath_c_dp.o : crmath.mod

