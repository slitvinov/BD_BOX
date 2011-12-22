#.cu.$(OBJEXT):
#	$(NVCC) -o $@ -c $< $(NVCCFLAGS)
.cu.$(OBJEXT):
	$(NVCC) -c -o $@ $< $(NVCCFLAGS) -Xcompiler "-I.. -DHAVE_CONFIG_H $(CFLAGS) $(OPENMP_CFLAGS)" -Xlinker "$(OPENMP_LDDFLAGS)" $(CUDPP_INCLUDE)
bin_PROGRAMS = bd_box
# set the include path found by configure
if CUDA
INCLUDES = $(all_includes) -I "$(NVCC_INCLUDE)"
else
INCLUDES = $(all_includes)
endif
bd_box_SOURCES = main.c \
                 data.c data.h \
                 input.c input.h \
                 output.c output.h \
                 err.c err.h \
                 rand_move.c rand_move.h \
                 trans.c trans.h \
                 math_help.h math_help.c \
                 main_loop.h main_loop.c \
                 myblas.h myblas.c \
                 restart.h restart.c \
                 potentials/bond.c potentials/bond.h \
                 potentials/angle.c potentials/angle.h \
                 potentials/angle_cos.c potentials/angle_cos.h \
                 potentials/dihe.c potentials/dihe.h \
                 potentials/dihe_angle.c potentials/dihe_angle.h \
                 potentials/electro.c potentials/electro.h \
                 potentials/electro_ext.c potentials/electro_ext.h \
                 potentials/LJ.c potentials/LJ.h \
                 potentials/pre_func.c \
                 potentials/bucket.c potentials/bucket.h \
                 potentials/sboundary.c potentials/sboundary.h \
                 diff_alg/cholesky_mpi.c diff_alg/cholesky_mpi.h \
                 diff_alg/newton.h \
                 diff_alg/diff_tensor.h \
                 potentials/calc_func.h \
                 cuda.c cuda.h \
                 potentials/sse_mathfun.h
if CUDA
bd_box_SOURCES += potentials/cuda/calc_func.cu
bd_box_SOURCES += potentials/cuda/sorter.cu
bd_box_SOURCES += potentials/cuda/sorter.cuh
bd_box_SOURCES += diff_alg/cuda/diff_tensor.cu
bd_box_SOURCES += diff_alg/cuda/diff_tensor.cuh
bd_box_SOURCES += diff_alg/cuda/diff_rand.cu
bd_box_SOURCES += diff_alg/cuda/diff_comp.cu
bd_box_SOURCES += diff_alg/cuda/diff_decl.cu
bd_box_SOURCES += diff_alg/cuda/diff_equation.cu
bd_box_SOURCES += diff_alg/cuda/newton.cu
bd_box_SOURCES += netdecl.cuh
else
if SSE
bd_box_SOURCES += potentials/calc_func_sse.c
bd_box_SOURCES += diff_alg/diff_tensor.c
bd_box_SOURCES += diff_alg/diff_rand.c
bd_box_SOURCES += diff_alg/diff_comp.c
bd_box_SOURCES += diff_alg/diff_decl.c
bd_box_SOURCES += diff_alg/diff_equation.c
bd_box_SOURCES += diff_alg/newton.c
else
bd_box_SOURCES += potentials/calc_func.c
bd_box_SOURCES += diff_alg/diff_tensor.c
bd_box_SOURCES += diff_alg/diff_rand.c
bd_box_SOURCES += diff_alg/diff_comp.c
bd_box_SOURCES += diff_alg/diff_decl.c
bd_box_SOURCES += diff_alg/diff_equation.c
bd_box_SOURCES += diff_alg/newton.c
endif
endif
bd_box_LDFLAGS = $(OPENMP_LDDFLAGS)
AM_CFLAGS = $(OPENMP_CFLAGS)
