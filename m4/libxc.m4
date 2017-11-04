AC_DEFUN([ACX_LIBXC], [
acx_libxc_ok=no

acx_libxc_save_CFLAGS="$CFLAGS"
acx_libxc_save_LIBS="$LIBS"

AC_ARG_WITH(libxc-prefix,
	[AC_HELP_STRING([--with-libxc-prefix], [directory where the libxc library is installed])])

LIBS="$acx_libxc_save_LIBS"
CFLAGS="$acx_libxc_save_CFLAGS"

if test x"$with_libxc_prefix" != x; then
  CFLAGS_LIBXC="-I$with_libxc_prefix/include"
  LIBS_LIBXC="-L$with_libxc_prefix/lib"
else
  LIBS_LIBXC=""
fi

CFLAGS="$CFLAGS_LIBXC $CFLAGS"
LIBS="$LIBS_LIBXC $LIBS $FCLIBS"

AC_CHECK_LIB(xc, xc_func_info_get_name, [acx_libxc_ok=yes], [])

if test x"$acx_libxc_ok" == xyes ; then
  LIBS_LIBXC="$LIBS_LIBXC -lxc"
  AC_DEFINE(HAVE_LIBXC, 1, [Define if you have a LIBXC library.])
fi

LIBS="$acx_libxc_save_LIBS"
CFLAGS="$acx_libxc_save_CFLAGS"

])dnl ACX_LIBXC
