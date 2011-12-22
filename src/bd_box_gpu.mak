!include <Win32.Mak>
.SUFFIXES: .exe .obj .c .cpp .cu

cflags = $(cflags) /I "$(CUDA_INC_PATH)"
conlibs = $(conlibs) /LIBPATH:"$(CUDA_LIB_PATH)" /LIBPATH:"$(NVSDKCOMPUTE_ROOT)\C\common\lib"
!IF "$(CPU)" == "i386"
conlibs = $(conlibs) cudart.lib cublas.lib cudpp32.lib
!ELSE
conlibs = $(conlibs) cudart.lib cublas.lib cudpp64.lib
!ENDIF
cvars = $(cvars) -DUSE_CUDA -DUSE_FLOAT
#for 64 also WIN32!!
cnvccvars = $(cvars) -DUSE_CUDPP -DWIN32 -DHAVE_ISINF -DHAVE_ISNAN
cnvccinc = -I"$(NVSDKCOMPUTE_ROOT)\C\common\inc\cudpp"
#change for your architecture
cnvccopt = -O3 -arch=sm_11
nvcc = "$(CUDA_BIN_PATH)\nvcc"
copenmp = 
!IFDEF OPENMP
copenmp=/openmp
!ENDIF

all: $(OUTDIR) $(OUTDIR)\bd_box.exe

$(OUTDIR) :
    if not exist "$(OUTDIR)/$(NULL)" mkdir $(OUTDIR)
	
.c{$(OUTDIR)}.obj:
    $(cc) $(cflags) $(cdebug) $(cvars) $(copenmp) /TC /Fo"$(OUTDIR)\\" /Fd"$(OUTDIR)\\" $<
	
{win\}.cpp{win\}.obj:
    $(cc) $(cflags) $(cdebug) $(cvars) $(copenmp) /TP /EHsc /Fo"win\\" /Fd"win\\" $<
CS = diff_alg\cuda
PS = potentials\cuda
{$(CS)\}.cu{$(CS)\}.obj:
    $(nvcc) $(cnvccopt) $(cnvccinc) -Xcompiler "$(cdebug) $(cnvccvars) $(copenmp)" -c $< -o $@ > nul

{$(PS)\}.cu{$(PS)\}.obj: 
    $(nvcc) $(cnvccopt) $(cnvccinc) -Xcompiler "$(cdebug) $(cnvccvars) $(copenmp)" -c $< -o $@ > nul
	
{potentials\}.c{potentials\}.obj:
    $(cc) $(cflags) $(cdebug) $(cvars) $(copenmp) /TC /Fo"potentials\\" /Fd"potentials\\" $<
{diff_alg\}.c{diff_alg\}.obj::
    $(cc) $(cflags) $(cdebug) $(cvars) $(copenmp) /TC /Fo"diff_alg\\" /Fd"diff_alg\\" $<
	
SS = $(CS)\diff_comp.obj $(CS)\diff_decl.obj $(CS)\diff_equation.obj \
 $(CS)\diff_rand.obj $(CS)\diff_tensor.obj $(CS)\newton.obj diff_alg\cholesky_mpi.obj \
 $(PS)\calc_func.obj $(PS)\sorter.obj

$(OUTDIR)\bd_box.exe: $(OUTDIR)\*.obj win\*.obj $(SS) potentials\*.obj
    $(link) $(conflags) -out:$(OUTDIR)\bd_box.exe $** $(conlibs)
	
clean:
        $(CLEANUP)
		del /Q potentials\cuda\*.obj
		del /Q potentials\*.obj
		del /Q diff_alg\cuda\*.obj
		del /Q diff_alg\*.obj
		del /Q win\*.obj
