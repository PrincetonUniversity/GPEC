# This is the umbrella makefile for stride_1.0.

# When first installed, we have to unpack the tarball (dislike)
SLU_COMMAND = cd SuperLU; ./remakeSuperLU.sh
ifneq ($(wildcard ./SuperLU/SuperLU_5.2.1/.),)
	SLU_COMMAND = cd SuperLU/SuperLU_5.2.1; make
endif

stride:
	$(SLU_COMMAND)
	cd SuperLU/SuperLU_5.2.1; make
	cd zvode; make
	cd zlange; make
	cd lsode; make
	cd vacuum; make
	cd equil; make
	cd dcon; make
	cp -f dcon/stride ../bin/

clean:
	rm -f *~
	cd SuperLU/SuperLU_5.2.1; make clean
	cd zvode; make clean
	cd zlange; make clean
	cd lsode; make clean
	cd vacuum; make clean
	cd equil; make clean
	cd dcon; make clean