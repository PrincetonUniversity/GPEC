#Global makefile for the DCON/IPEC suite of codes.

.IGNORE:

all_x: lsode_x equil_x orbit_x vacuum_x dcon_x match_x multi_x sum_x \
	xdraw_x ipec_x rundir_x

all: all_x

# executables_x:

xdraw_x:
	cd xdraw; make -f makefile_`uname`

lsode_x:
	cd lsode; make -f makefile_`uname`$(FORTRAN)

equil_x:
	cd equil; make -f makefile_`uname`$(FORTRAN)

orbit_x:
	cd orbit; make -f makefile_`uname`$(FORTRAN)

vacuum_x:
	cd vacuum; make -f makefile_`uname`$(FORTRAN)

dcon_x:
	cd dcon; make -f makefile_`uname`$(FORTRAN)

match_x:
	cd match; make -f makefile_`uname`$(FORTRAN)

multi_x:
	cd multi; make -f makefile_`uname`$(FORTRAN)

sum_x:
	cd sum; make -f makefile_`uname`$(FORTRAN)

ipec_x: 
	cd ipec; make -f makefile_`uname`$(FORTRAN)

rundir_x:
	mkdir -p rundir/`uname`$(FORTRAN)
	cp -f orbit/orbit rundir/`uname`$(FORTRAN)
	cp -f dcon/dcon rundir/`uname`$(FORTRAN)
	cp -f xdraw/xdraw rundir/`uname`$(FORTRAN)
	cp -f match/match rundir/`uname`$(FORTRAN)
	cp -f multi/multi rundir/`uname`$(FORTRAN)
	cp -f sum/sum rundir/`uname`$(FORTRAN)
	cp -f ipec/ipec rundir/`uname`$(FORTRAN)
	cp -f input/*.in* rundir/`uname`$(FORTRAN)
	cp -f input/mput rundir/`uname`$(FORTRAN)
	cp draw/drawdcon.in rundir/`uname`$(FORTRAN)
	cp draw/drawcrit.in rundir/`uname`$(FORTRAN)
	cp draw/drawsum.in rundir/`uname`$(FORTRAN)
	cp draw/draworbit.in rundir/`uname`$(FORTRAN)
	cp draw/drawprof.in rundir/`uname`$(FORTRAN)
	cp draw/draw2d.in rundir/`uname`$(FORTRAN)
	cp draw/drawahb.in rundir/`uname`$(FORTRAN)
	cp draw/drawgsei.in rundir/`uname`$(FORTRAN)
	cp draw/drawsol.in rundir/`uname`$(FORTRAN)
	cp draw/drawcon.in rundir/`uname`$(FORTRAN)
	cp draw/drawsurf.in rundir/`uname`$(FORTRAN)
	cp draw/drawdiffw.in rundir/`uname`$(FORTRAN)

clean:
	cd xdraw; make -f makefile_`uname` clean
	cd lsode; make -f makefile_`uname`$(FORTRAN) clean
	cd equil; make -f makefile_`uname`$(FORTRAN) clean
	cd orbit; make -f makefile_`uname`$(FORTRAN) clean
	cd vacuum; make -f makefile_`uname`$(FORTRAN) clean
	cd dcon; make -f makefile_`uname`$(FORTRAN) clean
	cd match; make -f makefile_`uname`$(FORTRAN) clean
	cd multi; make -f makefile_`uname`$(FORTRAN) clean
	cd sum; make -f makefile_`uname`$(FORTRAN) clean
	cd ipec; make -f makefile_`uname`$(FORTRAN) clean
	cd tex/orbit; rm -f *.dvi *.log
	cd tex/vacuum; rm -f *.dvi *.log
	cd tex/dcon/manuscript; rm -f *.dvi *.log
	cd tex/dcon/notebook; rm -f *.dvi *.log

realclean: clean
	rm -rf lib rundir
