AC_DEFUN([ACX_LAPACK], 
[
acx_lapack_ok=no

AC_ARG_WITH(lapack,
	[AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
	yes | "") ;;
	no) acx_lapack_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
	*) LAPACK_LIBS="-l$with_lapack" ;;
esac

acx_lapack_save_LIBS="$LIBS"

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
# Get fortran symbol name of dpotrf (a lapack routine).
	dpotrf=dpotrf_
	save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $LIBS"
	AC_MSG_CHECKING([for $dpotrf in $LAPACK_LIBS])
	AC_TRY_LINK_FUNC($dpotrf, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
	AC_MSG_RESULT($acx_lapack_ok)
	LIBS="$save_LIBS"
fi
# restore the libs variable
LIBS=$acx_lapack_save_LIBS
]) 
