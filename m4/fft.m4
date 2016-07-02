AC_DEFUN([ACX_FFTW3], [
AC_PREREQ(2.50)
acx_fftw3_ok=no

acx_fftw3_save_CFLAGS="$CFLAGS"
acx_fftw3_save_LIBS="$LIBS"

AC_ARG_WITH(fftw3-prefix,
	[AC_HELP_STRING([--with-fftw3-prefix], [directory where the ossp-fftw3 library is installed])])

if test x"$with_fftw3_prefix" != x; then
   CFLAGS_FFTW3="-I$with_fftw3_prefix/include"
   LIBS_FFTW3="-L$with_fftw3_prefix/lib -lfftw3"
else
   LIBS_FFTW3="-lfftw3"
fi

CFLAGS="$CFLAGS_FFTW3 $CFLAGS"
LIBS="$LIBS_FFTW3 $LIBS"

AC_CHECK_FUNC(fftw_plan_dft_1d, [acx_fftw3_ok=yes], [])

AC_MSG_CHECKING(for fftw3)
AC_MSG_RESULT($acx_fftw3_ok ( $CFLAGS_FFTW3 $LIBS_FFTW3 ))

if test x"$acx_fftw3_ok" == xyes ; then
  AC_DEFINE(HAVE_FFTW3, 1, [Define if you have a FFTW3 library.])
fi

LIBS="$acx_fftw3_save_LIBS"
CFLAGS="$acx_fftw3_save_CFLAGS"

])dnl ACX_FFTW3

dnl -----------------------------------------------------------------

AC_DEFUN([ACX_FFTW2], [
AC_PREREQ(2.50)
acx_fftw2_ok=no

acx_fftw2_save_CFLAGS="$CFLAGS"
acx_fftw2_save_LIBS="$LIBS"

AC_ARG_WITH(fftw2-prefix,
	[AC_HELP_STRING([--with-fftw2-prefix], [directory where the ossp-fftw2 library is installed])])

if test x"$with_fftw2_prefix" != x; then
   CFLAGS_FFTW2="-I$with_fftw2_prefix/include"
   LIBS_FFTW2="-L$with_fftw2_prefix/lib -lfftw"
else
   LIBS_FFTW2="-lfftw"
fi

CFLAGS="$CFLAGS_FFTW2 $CFLAGS"
LIBS="$LIBS_FFTW2 $LIBS"

AC_CHECK_FUNC(fftw_create_plan, [acx_fftw2_ok=yes], [])

AC_MSG_CHECKING(for fftw2)
AC_MSG_RESULT($acx_fftw2_ok ( $CFLAGS_FFTW2 $LIBS_FFTW2 ))

if test x"$acx_fftw2_ok" == xyes ; then
  AC_DEFINE(HAVE_FFTW2, 1, [Define if you have a FFTW2 library.])
fi

LIBS="$acx_fftw2_save_LIBS"
CFLAGS="$acx_fftw2_save_CFLAGS"

])dnl ACX_FFTW2
