!include <Win32.Mak>

copenmp = 
!IFDEF OPENMP
copenmp=/openmp
!ENDIF

all: $(OUTDIR) $(OUTDIR)\bd_box.exe

$(OUTDIR) :
    if not exist "$(OUTDIR)/$(NULL)" mkdir $(OUTDIR)
	
.c{$(OUTDIR)}.obj::
    $(cc) $(cflags) $(cdebug) $(cvars) $(copenmp) /TC /Fo"$(OUTDIR)\\" /Fd"$(OUTDIR)\\" $<
	
{win\}.cpp{win\}.obj::
    $(cc) $(cflags) $(cdebug) $(cvars) $(copenmp) /TP /EHsc /Fo"win\\" /Fd"win\\" $<

{diff_alg\}.c{diff_alg\}.obj::
    $(cc) $(cflags) $(cdebug) $(cvars) $(copenmp) /TC /Fo"diff_alg\\" /Fd"diff_alg\\" $<

{potentials\}.c{potentials\}.obj:
    $(cc) $(cflags) $(cdebug) $(cvars) $(copenmp) /TC /Fo"potentials\\" /Fd"potentials\\" $<
	
$(OUTDIR)\bd_box.exe: $(OUTDIR)\*.obj win\*.obj diff_alg\*.obj potentials\*.obj
    $(link) $(conflags) -out:$(OUTDIR)\bd_box.exe $** $(conlibs)
	
clean:
        $(CLEANUP)
		del /Q potentials\*.obj
		del /Q diff_alg\*.obj
		del /Q win\*.obj
