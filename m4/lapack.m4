dnl Available from the GNU Autoconf Macro Archive at:
dnl http://www.gnu.org/software/ac-archive/htmldoc/acx_lapack.html
dnl
AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

dnl We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
  acx_lapack_ok=noblas
fi

dnl Get fortran linker name of LAPACK function to check for.
dnl if not compiling with fortran, convert the names
m4_if(_AC_LANG, Fortran, [cheev=cheev], [AC_FC_FUNC(cheev)])

dnl Check if the library was given in the command line
if test $acx_lapack_ok = no; then
  AC_ARG_WITH(lapack, [AS_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
  case $with_lapack in
    yes | "") ;;
    no) acx_lapack_ok=disable ;;
    -* | */* | *.a | *.so | *.so.* | *.o) LIBS_LAPACK="$with_lapack" ;;
    *) LIBS_LAPACK="-l$with_lapack" ;;
  esac
fi

dnl Backup LIBS 
acx_lapack_save_LIBS="$LIBS"
LIBS="$LIBS_LAPACK $LIBS_BLAS $LIBS $FLIBS"

dnl First, check LIBS_LAPACK environment variable
if test $acx_lapack_ok = no; then
  AC_MSG_CHECKING([for $cheev in $LIBS_LAPACK])
  AC_LINK_IFELSE([AC_LANG_CALL([], [$cheev])], [acx_lapack_ok=yes], [])
  if test $acx_lapack_ok = no; then
    AC_MSG_RESULT([$acx_lapack_ok])
  else
    AC_MSG_RESULT([$acx_lapack_ok ($LIBS_LAPACK)])
  fi
fi

dnl Generic LAPACK library?
for lapack in mkl_lapack lapack lapack_rs6k acml; do
  if test $acx_lapack_ok = no; then
    AC_CHECK_LIB($lapack, $cheev,
      [acx_lapack_ok=yes; LIBS_LAPACK="$LIBS_LAPACK -l$lapack"], [], [$FLIBS])
  fi
done

AC_SUBST(LIBS_LAPACK)
LIBS="$acx_lapack_save_LIBS"

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
  AC_DEFINE(HAVE_LAPACK,1,[Defined if you have LAPACK library.])
  $1
else
  $2
fi
])dnl ACX_LAPACK
