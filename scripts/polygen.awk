#! /usr/bin/awk -f
# generate an initial configuration of the polyme
BEGIN {
# polymer
    # the number of beads
    N=60
    # distance between neghiboring beads
    dy=25
    # the hydrodynamic radius of the subunit
    sigma=12.25
    # the central charge
    Q=0.0
    # the hard-core (Lennard-Jones) radius (note that the value of R
    # should be doubled in the structure file, 2R)
    R=24.5
    # Lennard-Jones well depth
    eps=1.7766
    # mass
    mass=1

    # the equilibrium bond length
    r0=0.0
    # the maximum bond length
    rmax=36.75
    # the force constan
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
