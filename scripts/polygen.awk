BEGIN {

# polymer
    N=60
    dy=25
    sigma=12.25
    Q=0.0
    R=24.5
    eps=1.7766
    mass=1

# bond
    r0=24.5
    rmax=2.5e+07
    H=0.0986589

    x0=0.0
    y0=-N/2.0*dy
    z0=0.0

    for (id=1; id<N+1; id++) {
	printf("sub DNA %i %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n", 
	       id, x0, y0, z0, sigma, Q, R, eps, mass)
	y0+=dy
    }

    for (id=1; id<N; id++) {
	printf("bond %i %i %.2f %.2f %.2f\n",
	       id, id+1, r0, rmax, H);
    }
}
