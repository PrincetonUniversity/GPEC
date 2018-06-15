# This is the makefile for EQUIL.

#FFLAGS = -xHost -g -check all -openmp
#    ^----------Works for intel/2013.1.3
#FFLAGS = -qopt-report -xHost -g -check all -qopenmp -mkl -diag-disable 10397
#FFLAGS = -qopt-report -xHost -O3 -qopenmp -mkl -diag-disable 10397 #-pg #-ipo
FFLAGS = -xHost -O3 -qopenmp -mkl -diag-disable 10397 #-qopt-zmm-usage=high #-ipo

F90 = ifort $(FFLAGS)

CFLAGS = -O3 -DNDEBUG -DPRNTlevel=0 -DAdd_
CC = gcc $(CFLAGS)

HEADER = ../SuperLU/SuperLU_5.2.1/SRC

OBJS = \
	io.o \
	local.o \
	utils.o \
	spline.o \
	bicube.o \
	cspline.o \
	fft.o \
	c_fortran_zgssv.o \
	sparse.o \
	fspline.o \
	global.o \
	inverse.o \
	direct.o \
	read_eq.o \
	lar.o \
	gsec.o \
	sol.o \
	equil.o \
	equil_out.o

.f.o:
	$(F90) -c $*.f

.c.o:
	$(CC) -I$(HEADER) -c $<

libequil: $(OBJS)
	mkdir -p ../lib
	ar -r ../lib/libequil.a *.o

# dependencies

local.o: io.o
utils.o: local.o
spline.o: local.o
bicube.o: spline.o
cspline.o: local.o
fft.o: local.o
c_fortran_zgssv.o:
sparse.o: local.o c_fortran_zgssv.o
fspline.o: local.o spline.o cspline.o utils.o fft.o
global.o: bicube.o
inverse.o: global.o
direct.o: global.o
read_eq.o: inverse.o direct.o
lar.o: inverse.o
gsec.o: inverse.o
sol.o: direct.o
equil.o: read_eq.o gsec.o sol.o lar.o
equil_out.o: global.o

clean:
	rm -f *.o *.mod *~
