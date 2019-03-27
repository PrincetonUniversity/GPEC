#!/bin/bash
#THIS REGENERATES THE FORTRAN FILES FOR SuperLU
#http://people.sc.fsu.edu/~jburkardt/f_src/superlu/superlu.html
#Namely, it focuses on compiling and running: d_sample.f90
#RUN THIS SCRIPT:

rm *~
rm -rf SuperLU_5.2.1
tar -zxvf superlu_5.2.1.tar.gz
cp make.inc SuperLU_5.2.1/make.inc
sed -i 's@CURDIR@'"$PWD"'@' SuperLU_5.2.1/make.inc
cp Makefile SuperLU_5.2.1/Makefile
#Note: libblas.a NEEDN'T BE COPIED BY cp HERE,
#      (LIKE THE 'make' FILES ABOVE ARE)
#      BUT IS REFERENCED IN make.inc,
#      AND THEREFORE MUST BE KEPT!
#      NOTE IT WAS ORIGINALLY COPIED FROM pppl.

cd SuperLU_5.2.1
make

#Then to run something in Fortran...
cd FORTRAN
cp ../../F90_Code/d_sample.f90 d_sample.f90
cp ../../F90_Code/Makefile Makefile
make clean
make
#./d_sample

#If desired to work with SuperLU's code directly,
#and not external wrappers, then run:
#./df77exm < ../EXAMPLE/g20.rua
