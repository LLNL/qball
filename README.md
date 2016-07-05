# Qball

## Installing

To compile Qbox:

1. Go to trunk/src

2. Create an architecture include file, e.g. bgq_ctf_essl.mk.  This file should
define the make variables for compiling and linking (CXX, CXXFLAGS, LDFLAGS
etc.) and should contain the locations of all libraries and include files.

3. Qbox requires ScaLAPACK 2 and dependent libraries (blas and lapack
or vendor-equivalents, e.g. ESSL), FFTW 2.x or equivalent, and
(optionally) Xerces XML. You also need the makedepend utility (in
Debian comes in the xutils-dev package).

4. Build the code using the ARCH variable to point to this include file, e.g.:
    
   ```
   make ARCH=bgq_ctf_essl
   ```
   
   If the ARCH variable is not specified, a default value will be guessed
   from the hostname or uname commands listed in Makefile.arch.  

5. Object files are stored in objs-$(ARCH) subdirectories, so one can
compile multiple versions of the code simultaneously in the same
directory.

Contact Erik Draeger (draeger1@llnl.gov) or Xavier Andrade
(xavier@llnl.gov) with any questions or problems.

## Running

To run Qbox, one needs an input file (.i), a coordinate file (.sys)
and pseudopotential file(s) (.xml).  Input examples can be found in
the examples/ directory.

The input file can be specified either as an argument or as stdin to Qbox, e.g.

    srun -n 16384 ./qb-bgq_ctf_essl gold.N992.i > gold.N992.out

    srun -n 64 ./qb-linux < test.i > test.out

## Release

Qball is licensed under the terms of the [GPL v3 License](/COPYING).

``LLNL-CODE-635376``
