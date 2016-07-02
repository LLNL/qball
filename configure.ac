#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT(qball, alpha, xavier@llnl.gov)
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/FoldInWsCmd.h])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CXX(mpic++)
AC_PROG_CC(mpicc)
AC_PROG_CC_C99
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Checks for OpenMP
AC_LANG_PUSH(C++)
AX_OPENMP
CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
AC_LANG_POP(C++)

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FSEEKO
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday memset mkdir pow sqrt strspn uname])

#LIBRARIES

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR FFT LIBRARIES
================================================================================])

found_fft=no

#FFTW3
if test x$found_fft != xyes ; then
  ACX_FFTW3
fi

if test x$acx_fftw3_ok == xyes ; then
  found_fft=yes
  CFLAGS="$CFLAGS $CFLAGS_FFTW3"
  LIBS="$LIBS $LIBS_FFTW3"
fi

#FFTW2
if test x$found_fft != xyes ; then
   ACX_FFTW2
fi

if test x$acx_fftw2_ok == xyes ; then
  found_fft=yes
  CFLAGS="$CFLAGS $CFLAGS_FFTW2"
  LIBS="$LIBS $LIBS_FFTW2"
fi

if test x$found_fft != xyes ; then
  AC_MSG_WARN([

   Cannot find an FFT library, using the slower internal FFT routines.

  ])
  AC_DEFINE(FFT_NOLIB, 1, [no FFT library])
fi

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR MPI
================================================================================])

acx_mpi_ok=yes
AC_LANG_PUSH(C++)
AC_CHECK_HEADERS([mpi.h], [], [acx_mpi_ok=no])
AC_CHECK_FUNCS([MPI_Init], [], [acx_mpi_ok=no])
AC_LANG_POP(C++)
AC_MSG_CHECKING(whether we can use MPI)
AC_MSG_RESULT($acx_mpi_ok)

if test $acx_mpi_ok == yes; then
   AC_DEFINE(USE_MPI, 1, [enable MPI support])
else
   AC_MSG_ERROR([
   
   Cannot find the required MPI library. If you have an MPI Implementation
   installed you might need to set CC=mpicc and CXX=mpic++.
   ])
fi

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR THE OPTIONAL XERCES-C LIBRARY (USED TO PARSE XML)
================================================================================])
AX_LIB_XERCES
if test x$HAVE_XERCES == xyes ; then
   CFLAGS="$CFLAGS $XERCES_CFLAGS"
   LDFLAGS="$LDFLAGS $XERCES_LDFLAGS"
   LIBS="$LIBS $XERCES_LIBS"
   AC_DEFINE(USE_XERCES, 1, [have the Xerces-C library])
fi

AC_DEFINE(ARCH, "autotools", [architecture])

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR REQUIRED LINEAR ALGEBRA LIBRARIES (BLAS AND LAPACK)
================================================================================])
dnl We need to check for fortran, to link with blas and other libraries
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS

dnl check for BLAS
ACX_BLAS([], AC_MSG_ERROR([could not find required BLAS library]))

dnl check for LAPACK
ACX_LAPACK([], AC_MSG_ERROR([could not find required LAPACK library]))

AC_MSG_NOTICE([
================================================================================
  CHECKING FOR OPTIONAL PARALLEL LINEAR ALGEBRA LIBRARIES (BLACS AND SCALAPACK)
================================================================================])

ACX_BLACS
ACX_SCALAPACK
ACX_SCALAPACK2

LIBS="$LIBS_SCALAPACK $LIBS_BLACS $LIBS_LAPACK $LIBS_BLAS $LIBS $FCLIBS"

AC_MSG_NOTICE([
================================================================================
  CONFIGURATION DONE, GENERATING MAKEFILES
================================================================================])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 testsuite/Makefile])

AC_OUTPUT
