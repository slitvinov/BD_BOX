#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.64])
AC_INIT( BD_BOX, 1.0.2, mdlugosz@icm.edu.pl, bd_box.tgz, [http://bionano.icm.edu.pl/software/bd_box] )
AC_CONFIG_SRCDIR([src/input.h])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE(bd_box, 1.0.2)

# Checks for programs.
AC_PROG_CC

#use mpi
AC_ARG_ENABLE([mpi],
              [AS_HELP_STRING([--enable-mpi], [use mpi (default is no)])],
              [use_mpi=$enable_mpi],
              [use_mpi=no])

if test "$use_mpi" == yes; then
    AC_DEFINE([USE_MPI],[1],[Use MPI])
else
    AC_OPENMP
fi

#use sse
AC_ARG_ENABLE([sse],
              [AS_HELP_STRING([--enable-sse], [use sse (default is no)])],
              [use_sse=$enable_sse],
              [use_sse=no])

AM_CONDITIONAL([SSE], [test x$use_sse == xyes])
if test x$use_sse = xyes; then
    AC_DEFINE([USE_SSE],[1],[Use vectorization])
fi

#single precision
AC_ARG_ENABLE([float],
              [AS_HELP_STRING([--enable-float], [use single-precision (default is double-precision)])],
              [use_float=$enable_float],
              [use_float=no])
if test x$use_float = xyes; then
    AC_DEFINE([USE_FLOAT],[1],[Use floats])
fi


#check for ScaLAPACK pdpotrf
sinclude([acx_scalapack.m4])
ACX_SCALAPACK
if test $acx_scalapack_ok == "yes"; then
    LIBS+=" $SCALAPACK_LIBS"
    AC_DEFINE([USE_SCALAPACK],[1],[Use SCALAPACK])
fi

#check for GSL
sinclude([acx_gsl.m4])
ACX_GSL
if test $acx_gsl_ok == "yes"; then
    LIBS+=" $GSL_LIBS"
    AC_DEFINE([USE_GSL],[1],[Use GSL])
fi

#check for CUDA
AC_ARG_VAR([NVCCFLAGS], [Flags for CUDA compiler])
AC_ARG_VAR([NVCCLIB], [Flags for CUDA compiler lib directory lib or lib64])
sinclude([acx_cuda.m4])
ACX_MAGMA
if test x$acx_magma_ok = xyes; then
    LIBS+=" $MAGMA_LIBS"
    AC_DEFINE([USE_MAGMA],[1],[Use magma])
fi
ACX_CUDPP
if test x$acx_cudpp_ok = xyes; then
    LIBS+=" -L $CUDPP_LIBS -lcudpp"
    AC_SUBST([CUDPP_INCLUDE],["-I $CUDPP_LIBS"])
    AC_DEFINE([USE_CUDPP],[1],[Use cudpp])
fi

ACX_CUDA
if test -z "$NVCCLIB"; then
    NVCCLIB=lib64
fi
if test x$acx_cuda_ok = xyes; then
    LIBS+=" -L $CUDA_LIBS/$NVCCLIB -lcuda -lcudart -lcublas"
    AC_SUBST([NVCC],[$CUDA_LIBS/bin/nvcc])
    AC_SUBST([NVCC_INCLUDE],[$CUDA_LIBS/include])
    AC_SUBST([NVCCFLAGS])
    AC_DEFINE([USE_CUDA],[1],[Use CUDA])
fi

AM_CONDITIONAL([CUDA], [test x$acx_cuda_ok = xyes])

#check for LAPACK dpotrf
sinclude([acx_lapack.m4])
ACX_LAPACK
if test x$acx_lapack_ok = xyes; then
    LIBS+=" $LAPACK_LIBS"
    AC_DEFINE([USE_LAPACK],[1],[Use LAPACK])
fi

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])

# Checks for header files.
AC_CHECK_HEADERS([ limits.h malloc.h memory.h stdlib.h string.h float.h ])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CHECK_FUNCS([floor memset sqrt strchr drand48_r pow strrchr erfc isnan isinf])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
