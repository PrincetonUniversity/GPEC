#Global makefile for the GPEC suite of codes.

.IGNORE:

all_x: lsode_x equil_x orbit_x vacuum_x dcon_x match_x multi_x sum_x \
	xdraw_x coil_x ipec_x pent_x pentrc_x rundir_x

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

coil_x:
	cd coil; make -f makefile_`uname`$(FORTRAN)

ipec_x:  
	cd ipec; make -f makefile_`uname`$(FORTRAN)

pent_x:  
	cd pent; make -f makefile_`uname`$(FORTRAN)

pentrc_x:  
	cd pentrc; make -f makefile_`uname`$(FORTRAN)

rundir_x:
	mkdir -p rundir/`uname`$(FORTRAN)
	cp -f orbit/orbit rundir/`uname`$(FORTRAN)
	cp -f dcon/dcon rundir/`uname`$(FORTRAN)
	cp -f xdraw/xdraw rundir/`uname`$(FORTRAN)
	cp -f match/match rundir/`uname`$(FORTRAN)
	cp -f multi/multi rundir/`uname`$(FORTRAN)
	cp -f sum/sum rundir/`uname`$(FORTRAN)
	cp -f ipec/ipec rundir/`uname`$(FORTRAN)
	cp -f pent/pent rundir/`uname`$(FORTRAN)
	cp -f pent/*.dat rundir/`uname`$(FORTRAN)
	cp -f pentrc/pentrc rundir/`uname`$(FORTRAN)
	cp -f input/*.in* rundir/`uname`$(FORTRAN)
	cp -f input/mput rundir/`uname`$(FORTRAN)
	cp -f draw/drawdcon.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawcrit.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawsum.in rundir/`uname`$(FORTRAN)
	cp -f draw/draworbit.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawprof.in rundir/`uname`$(FORTRAN)
	cp -f draw/draw2d.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawahb.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawgsei.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawsol.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawcon.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawsurf.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawpmodb.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawpmodb_2d.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawxbnormal.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawvbnormal.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawxbtangent.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawxbtangent_2d.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawpflux_re_2d.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawpflux_im_2d.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawbnormal_spectrum.in rundir/`uname`$(FORTRAN)
	cp -f draw/drawvbnormal_spectrum.in rundir/`uname`$(FORTRAN)
	cp -f coil/*.dat rundir/`uname`$(FORTRAN)
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
	cd coil; make -f makefile_`uname`$(FORTRAN) clean
	cd ipec; make -f makefile_`uname`$(FORTRAN) clean
	cd pent; make -f makefile_`uname`$(FORTRAN) clean
	cd pentrc; make -f makefile_`uname`$(FORTRAN) clean
	cd tex/orbit; rm -f *.dvi *.log
	cd tex/vacuum; rm -f *.dvi *.log
	cd tex/dcon/manuscript; rm -f *.dvi *.log
	cd tex/dcon/notebook; rm -f *.dvi *.log

realclean: clean
	rm -rf lib rundir
