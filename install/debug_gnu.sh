#!/bin/sh

export FC=gfortran
export CC=gcc
export CXX=g++
export FFLAGS="-g -O0 -fbacktrace -fdefault-real-8 -fdefault-double-8 \
    -fallow-argument-mismatch -lnetcdf -lnetcdff \
    -finit-real=nan"
export NETCDF_FORTRAN_HOME=/usr/include
export LAPACK_HOME=/usr/include
OMPFLAG="-fopenmp" make NETCDFLIBS="-lnetcdf -lnetcdff"
