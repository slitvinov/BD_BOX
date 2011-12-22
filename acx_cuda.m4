AC_DEFUN([ACX_CUDA], 
[
acx_cuda_ok=no
AC_ARG_WITH(cuda,
    [AC_HELP_STRING([--with-cuda=<path>], [use cudatoolkit path])])
case $with_cuda in
	yes) acx_cuda_ok=yes ;;
	no | "") acx_cuda_ok=no ;;
	-* | *.a | *.so | *.so.* | *.o) AC_MSG_ERROR([Apply PATH to cudatoolkit]) ;;
	*) CUDA_LIBS="$with_cuda"; acx_cuda_ok=yes ;;
esac

]) 

AC_DEFUN([ACX_MAGMA],
[
acx_magma_ok=no
AC_ARG_WITH(magma,
	[AC_HELP_STRING([--with-magma=<lib>], [use magma library <lib>])])
case $with_magma in
	yes) acx_magma_ok=yes ;;
	no | "") acx_magma_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) MAGMA_LIBS="$with_magma"; acx_magma_ok=yes ;;
	*) MAGMA_LIBS="-l$with_magma"; acx_magma_ok=yes ;;
esac

]) 

AC_DEFUN([ACX_CUDPP],
[
acx_cudpp_ok=no
AC_ARG_WITH(cudpp,
	[AC_HELP_STRING([--with-cudpp=<path>], [use cudpp path, containing libcudpp.a and cudpp.h])])
case $with_cudpp in
	yes) acx_cudpp_ok=yes ;;
	no | "") acx_cudpp_ok=disable ;;
	-* | *.a | *.so | *.so.* | *.o) AC_MSG_ERROR([Apply PATH to cudpp ]) ;;
	*) CUDPP_LIBS="$with_cudpp"; acx_cudpp_ok=yes ;;
esac

]) 


