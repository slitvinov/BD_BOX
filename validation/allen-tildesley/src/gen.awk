#/usr/bin/awk -f

BEGIN {
    Nx=4
    Ny=4
    Nz=2
    N=Nx*Ny*Nz
    dx=dy=dz=1
    print N

    for (i=0; i<Nx; i++) {
	for (j=0; j<Ny; j++) {
	    for (k=0; k<Nz; k++) {
		print i*dx, j*dy, k*dz
	    }
	}
    }
}