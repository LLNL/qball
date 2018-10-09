# Qball

## Installing

To compile Qbox:

1. Generate the configure script (you need autoconf and automake). 

  ```
  autoreconf -i
  ```

2. Now you need to determine how to run the configure script. Since Qball depends on some non-standard libraries you might need to set some environment variables and to add some flags to tell qball where to find those libraries.

  | Environment variable | Description             | Note                  |
  |----------------------|-------------------------|-----------------------|
  | CC                   | C compiler              | Default is mpicc      |
  | CXX                  | C++ compiler            | Default is mpic++     |
  | FC                   | Fortran compiler        | Default is mpif90. Used only to detect Fortran libraries.|
  | CFLAGS               | C compiler flags        |                       |
  | CXXFLAGS             | C++ compiler flags      |                       |
  | FCFLAGS              | Fortran compiler flags  | Used only to detect Fortran libraries.|
  | LDFLAGS              | Flags to add to the linker |                    |
  | LIBS                 | Extra libs add to linking  |                    |
  | LIBS_BLAS            | Compilation flags to add the blas library |     |


  | Flag                  | Value                          | Note                  |
  |-----------------------|--------------------------------|-----------------------|
  | --prefix=             | installation directory         | default is /usr/local |
  | --with-fftw3-prefix=  | path where fftw3 is installed  |                       |
  | --with-fftw2-prefix=  | path where fftw2 is installed  |                       |
  | --with-essl-prefix=   | path where the IBM ESSL library is installed |         |
  | --with-blas=          | path where the Blas library file is located  | you can also use LIBS_BLAS |
  | --with-lapack=        | path where the lapack library file is located |        |
  | --with-blacs=         | path where the blacs library file is located | you can also pass the location of scalapack |
  | --with-scalapack=     | path where the scalapack file is located |             |

  For example, for a Blue Gene/Q system, you configure script might look something like this:

  ```
  QBALLPREFIX=/usr/local/
  QBALLDEPS=$QBALLPREFIX/dependencies/
  export CC=mpixlc_r
  export CXX=mpixlcxx_r
  export FC=mpixlf95_r
  export LIBS_BLAS="-L/usr/local/tools/essl/5.1/lib/ -lesslsmpbg"
  export LIBS_BLACS="-lmpiblacs"
  export LDFLAGS="-qsmp=omp"
  export CFLAGS="-O3 -qsmp=omp -qarch=qp -qtune=qp"
  export CXXFLAGS="$CFLAGS -qlanglvl=extended0x"
  export FCFLAGS=$CFLAGS" -qxlf90=autodealloc -qessl -I$HOME/$xarch/fftw-3.3.4/include"
  ./configure --with-essl-prefix=/usr/local/tools/essl/5.1/ \
    --with-lapack=$QBALLDEPS/lib/liblapack.a --with-blacs=$QBALLDEPS/lib/libscalapack.a --prefix=$QBALLPREFIX
  ```

3. Run the configure script with the necessary flags:
  
  ```
  ./configure --prefix=... 
  ```

4. Now we are ready to build the code:
    
   ```
   make
   make install
   ```

Contact Erik Draeger (draeger1@llnl.gov) or Xavier Andrade
(xavier@llnl.gov) with any questions or problems.

## Running

To run Qball, one needs an input file (.i), a coordinate file (.sys)
and pseudopotential file(s) (.xml).  Input examples can be found in
the examples/ directory.

The input file can be specified either as an argument or as stdin to qball, e.g.

    srun -n 16384 qball gold.N992.i > gold.N992.out

    srun -n 64 qball < test.i > test.out

## Release

Qball is licensed under the terms of the [GPL v3 License](/COPYING).

``LLNL-CODE-635376``
