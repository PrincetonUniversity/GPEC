
F90 = $(FC) $(FFLAGS)

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

# targets

all : superlu libequil

# need this for STRIDE's new sparse.f
superlu:
	$(SLU_COMMAND)

libequil: $(OBJS)
	mkdir -p ../lib/
	ar -r ../lib/libequil.a *.o

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
global.o: bicube.o cspline.o
inverse.o: global.o utils.o
direct.o: global.o
read_eq.o: inverse.o direct.o
lar.o: inverse.o
gsec.o: inverse.o
sol.o: direct.o
equil.o: read_eq.o lar.o gsec.o sol.o
equil_out.o: global.o
