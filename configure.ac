#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
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

# Checks for libraries.
AC_DEFINE(USE_MPI, 1, [enable MPI support])
AC_DEFINE(FFT_NOLIB, 1, [no FFT library])

# Checks for header files.
AC_CHECK_HEADERS([stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FSEEKO
AC_FUNC_STRTOD
AC_CHECK_FUNCS([gettimeofday memset mkdir pow sqrt strspn uname])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 testsuite/Makefile])

AC_OUTPUT
