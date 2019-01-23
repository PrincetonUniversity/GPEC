
FFLAGS = -qopenmp
EXTRA_FFLAGS = -xHost -mkl -diag-disable 10397 #-03  #-qopt-zmm-usage=high #-ipo

F90 = $(FC) $(FFLAGS) $(EXTRA_FFLAGS)

CFLAGS = -O3 -DNDEBUG -DPRNTlevel=0 -DAdd_
CC = gcc $(CFLAGS)
HEADER = ../superlu/SuperLU_5.2.1/SRC

.f.o:
	$(F90) -c $*.f

.c.o:
	$(CC) -I$(HEADER) -c $<

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

libequil: $(OBJS)
	mkdir -p ../lib/
	ar -r ../lib/libequil.a *.o

# targets

clean:
	rm -f *.o *.mod

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
equil.o: read_eq.o lar.o gsec.o sol.o
equil_out.o: global.o
