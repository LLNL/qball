AC_DEFUN([ACX_MASSV], [
AC_PREREQ(2.50)
acx_massv_ok=no

acx_massv_save_CFLAGS="$CFLAGS"
acx_massv_save_LIBS="$LIBS"

AC_ARG_WITH(massv-prefix,
	[AC_HELP_STRING([--with-massv-prefix], [directory where the massv library is installed])])

LIBS="$acx_massv_save_LIBS"
CFLAGS="$acx_massv_save_CFLAGS"

if test x"$with_massv_prefix" != x; then
  CFLAGS_MASSV="-I$with_massv_prefix/include"
  LIBS_MASSV="-L$with_massv_prefix/lib"
else
  LIBS_MASSV=""
fi

CFLAGS="$CFLAGS_MASSV $CFLAGS"
LIBS="$LIBS_MASSV $LIBS $FCLIBS"

AC_CHECK_LIB(massv, vsincos, [acx_massv_ok=yes], [])

if test x"$acx_massv_ok" == xyes ; then
  LIBS_MASSV="$LIBS_MASSV -lmassv"
  AC_DEFINE(HAVE_MASSV, 1, [Define if you have a MASSV library.])
fi

LIBS="$acx_massv_save_LIBS"
CFLAGS="$acx_massv_save_CFLAGS"

])dnl ACX_MASSV
