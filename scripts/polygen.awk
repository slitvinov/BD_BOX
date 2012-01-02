#! /usr/bin/awk -f
# generate an initial configuration of the polyme
BEGIN {
# polymer
    # the number of beads
    if (length(N)==0) {
	N=60
    }
    # the hard-core (Lennard-Jones) radius (note that the value of R
    # should be doubled in the structure file, 2R)
    # this is a diameter
    R=24.5
    # initial istance between neghiboring beads
    dy=1.1*R
    # the central charge
    Q=0.0
    # Lennard-Jones well depth
    eps=3.7766
    # mass
    mass=1
    # the equilibrium bond length
    r0=0.0
    # the maximum bond length
    rmax=1.5*R
    # the force constan
    H=0.0986589
    # the hydrodynamic radius of the subunit
    sigma=12.25


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
