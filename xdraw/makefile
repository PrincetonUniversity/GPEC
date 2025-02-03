# This is the makefile for XDRAW on Linux computers.
include ../install/DEFAULTS.inc

CFLAGS = -DUNIX -O2 $(X11_INC_DIR)
CC += $(CFLAGS)

.c.o:
	$(CC) -c $*.c

# Linux and OSX standards
LIBDIR= \
	-L/usr/X11R6/lib64 \
	-L/opt/local/include \
	-L/opt/homebrew/include \
	$(X11_LIB_DIR)

LIBS = \
	$(LIBDIR) \
	-lm \
	-lX11

OBJS1 = \
	binread.o \
	event.o \
	gradient.o \
	limits.o \
	ps.o \
	redraw.o \
	setcolor.o \
	spline.o \
	spline1.o \
	xcontour.o \
	xdraw.o \
	xedit.o \
	xinit.o \
	xtools.o

OBJS2 = \
	ps.o \
	pspack.o

OBJS3 = tell.o

all: xdraw pspack

xdraw: $(OBJS1)
	$(CC) -o xdraw $(OBJS1) $(LIBS)
	mkdir -p ../bin
	cp -f xdraw ../bin

pspack: $(OBJS2)
	$(CC) -o pspack $(OBJS2)

tell: $(OBJS3)
	$(CC) -o tell $(OBJS3)


# dependencies

binread.o: curves.h gendefs.h xinit.h
event.o: curves.h device.h gendefs.h glx.h setcolor.h xcontour.h \
	xdraw.h xedit.h xinit.h xtools.h
limits.o: curves.h gendefs.h xinit.h xtools.h
ps.o: gendefs.h ps.h
pspack.o: pspack.h
redraw.o: curves.h gendefs.h ps.h setcolor.h spline.h \
	xcontour.h xdraw.h xedit.h xinit.h xtools.h
setcolor.o: gendefs.h ps.h setcolor.h
spline.o: gendefs.h spline.h
xcontour.o: curves.h gendefs.h ps.h setcolor.h xcontour.h xdraw.h \
	xtools.h
xdraw.o: curves.h gendefs.h ps.h setcolor.h wdraw.h xcontour.h \
	xdraw.h xinit.h xtools.h
xedit.o: curves.h gendefs.h xinit.h xtools.h
xinit.o: curves.h gendefs.h ps.h xcontour.h xdraw.h xedit.h xinit.h \
	xtools.h
xtools.o: curves.h gendefs.h ps.h setcolor.h xcontour.h xdraw.h xedit.h \
	xinit.h xtools.h

tell.o: tell.c

clean:
	rm -f *.o xdraw pspack
