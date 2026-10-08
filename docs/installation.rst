.. _installation:

Installation
============

This page describes how to build GPEC from source using the Makefiles in
the install directory. The build system supports multiple compilers and
platforms and can automatically build missing dependencies when needed.

Overview
--------

The top-level build entry point is the Makefile in the install directory.
It compiles all major modules and stages executables into the bin directory.

Key behaviors:

* Auto-detects math and netCDF libraries (MKL, OpenBLAS, LAPACK, netCDF).
* Builds missing dependencies into deps (OpenBLAS, HDF5, netCDF-C, netCDF-Fortran).
* Uses git metadata to generate version files for some modules.

Requirements
------------

* Fortran compiler (ifort/ifx, gfortran, or another modern OpenMP-capable compiler)
* C compiler (for xdraw and some dependencies)
* git (required for version metadata and to fetch dependency submodules)
* BLAS/LAPACK and netCDF (auto-built if not found)

Quick start (local build)
-------------------------

From the install directory, run::

	make

This builds all main executables and libraries. Output executables appear in
bin and also inside each module subdirectory.

To build a single module, use::

	make dcon
	make gpec

To print the detected configuration, run::

	make v

Important make targets
----------------------

* all (default): build all main executables and libraries
* allf: same as all but Fortran-only modules
* clean: remove built objects and module files
* clear: remove example outputs in docs/examples
* realclean: deep clean including dependencies built into deps
* <module>: build a specific module (e.g. dcon, gpec, stride)

Build options and environment variables
---------------------------------------

You can override compilers and library locations without editing files by
setting environment variables or make variables.

Common overrides:

* FC: Fortran compiler (e.g. ifort, ifx, gfortran)
* CC: C compiler (for xdraw and dependencies)
* FFLAGS: Fortran flags (optimization, debugging, OpenMP)
* CFLAGS: C compiler flags
* OMPFLAG: OpenMP flag (e.g. -qopenmp or -fopenmp)

Math libraries (BLAS/LAPACK):

* MKLROOT: use Intel MKL
* OPENBLASHOME: use OpenBLAS
* LAPACKHOME or LAPACK_HOME: use system LAPACK/BLAS
* ACML_HOME: use ACML (legacy)

netCDF libraries:

* NETCDF_FORTRAN_HOME: preferred when C and Fortran installs are split
* NETCDFHOME or NETCDF_DIR: netCDF installation prefix
* NETCDFINC: netCDF include path
* NETCDF_C_HOME or NETCDF_C_DIR: netCDF-C location if separate

Dependency install location:

* DEPSINSTALLDIR: where auto-built deps are installed (default: deps)

Examples:

Build with debug flags::

	make FFLAGS="-g -traceback -check all -recursive"

Build with modern gfortran (avoid argument mismatch errors)::

	make FFLAGS="-fallow-argument-mismatch"

Automatic dependency builds
---------------------------

If the build cannot find required math or netCDF libraries, it will build
dependencies into deps. These include:

* OpenBLAS
* HDF5
* netCDF-C
* netCDF-Fortran

To build all dependencies explicitly, run::

	make alldeps

Platform notes
--------------

Linux clusters (module environments)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

On supported clusters, it is recommended to use the public GPEC module
which sets compilers and library paths correctly. Example::

	module load gpec

To use the development module, load gpec/0.0 if provided. After loading a
module, verify configuration with::

	make v

Generic Linux (gfortran)
^^^^^^^^^^^^^^^^^^^^^^^^

Set paths manually if modules are unavailable::

	export FC=gfortran
	export LAPACKHOME=/usr/lib/x86_64-linux-gnu/
	export NETCDFHOME=/usr/lib/x86_64-linux-gnu/
	export NETCDFINC=/usr/include

Then build with::

	make FFLAGS=""

macOS (Intel)
^^^^^^^^^^^^^

Install gfortran, netCDF, and LAPACK (via MacPorts or Homebrew). Typical
settings are::

	export FC=/usr/local/bin/gfortran
	export NETCDFHOME=/opt/local
	export LAPACKHOME=/usr/lib

Then build with::

	make FFLAGS="-fallow-argument-mismatch"

macOS (Apple Silicon)
^^^^^^^^^^^^^^^^^^^^^

If you do not have the Xcode command line tools installed, run::

    xcode-select --install
    
Use Homebrew for netCDF and the Accelerate framework for BLAS/LAPACK. Example::

	export CC=/opt/homebrew/bin/gcc-14
	export FC=/opt/homebrew/bin/gfortran-14
	export NETCDFHOME=/opt/homebrew
    export X11_HOME=/opt/homebrew
	export LAPACKHOME=/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/

Then build with::

	make FFLAGS="-fallow-argument-mismatch"


X11 (xdraw)
^^^^^^^^^^^

If xdraw fails to link on macOS, ensure XQuartz is installed and add include
and library paths in xdraw/makefile, for example::

    brew install libx11
    export X11_HOME=/opt/homebrew

Troubleshooting
---------------

* Run make v to verify detected compiler and library paths.
* If netCDF or BLAS/LAPACK are not found, set NETCDFHOME and LAPACKHOME
  (or MKLROOT/OPENBLASHOME) explicitly.
* Ensure git is available so version files can be generated.

If your platform is not covered above, provide your OS, compiler, and library
installation paths to a maintainer and update DEFAULTS.inc as needed.
