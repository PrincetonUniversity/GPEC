# This is the umbrella makefile for stride_1.0.

stride:
	cd SuperLU; ./remakeSuperLU.sh
	cd zvode; make
	cd zlange; make
	cd lsode; make
	cd vacuum; make
	cd equil; make
	cd dcon; make

clean:
	rm -f *~
	cd SuperLU/SuperLU_5.2.1; make clean
	cd zvode; make clean
	cd zlange; make clean
	cd lsode; make clean
	cd vacuum; make clean
	cd equil; make clean
	cd dcon; make clean