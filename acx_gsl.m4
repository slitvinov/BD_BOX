AC_DEFUN([ACX_GSL], 
[
acx_gsl_ok=no

AC_ARG_WITH(gsl,
	[AC_HELP_STRING([--with-gsl=<lib>], [use gsl library <lib>])])
case $with_gsl in
	yes) acx_gsl_ok=yes ;;
	no | "") acx_gsl_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) GSL_LIBS="$with_gsl"; acx_gsl_ok=yes ;;
	*) GSL_LIBS="-l$with_gsl"; acx_gsl_ok=yes ;;
esac

]) 
