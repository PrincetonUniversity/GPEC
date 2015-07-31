#Global makefile for the GPEC suite of codes.
 
.IGNORE:

all: lsode_x equil_x orbit_x vacuum_x dcon_x match_x multi_x sum_x \
	xdraw_x coil_x ipec_x pentrc_x bin_x

# executables_x:

xdraw_x:
	cd xdraw; make -f makefile

lsode_x:
	cd lsode; make -f makefile

equil_x:
	cd equil; make -f makefile

orbit_x:
	cd orbit; make -f makefile

vacuum_x:
	cd vacuum; make -f makefile

dcon_x:
	cd dcon; make -f makefile

match_x:
	cd match; make -f makefile

multi_x:
	cd multi; make -f makefile

sum_x:
	cd sum; make -f makefile

coil_x:
	cd coil; make -f makefile

ipec_x:  
	cd ipec; make -f makefile

pentrc_x:  
	cd pentrc; make -f makefile

bin_x:
	mkdir -p bin
	cp -f orbit/orbit bin
	cp -f dcon/dcon bin
	cp -f xdraw/xdraw bin
	cp -f match/match bin
	cp -f multi/multi bin
	cp -f sum/sum bin
	cp -f ipec/ipec bin
	cp -f pentrc/pentrc bin
	cp -f input/mput bin

clean:
	cd xdraw; make -f makefile clean
	cd lsode; make -f makefile clean
	cd equil; make -f makefile clean
	cd orbit; make -f makefile clean
	cd vacuum; make -f makefile clean
	cd dcon; make -f makefile clean
	cd match; make -f makefile clean
	cd multi; make -f makefile clean
	cd sum; make -f makefile clean
	cd coil; make -f makefile clean
	cd ipec; make -f makefile clean
	cd pentrc; make -f makefile clean

realclean: clean
	rm -rf lib bin
