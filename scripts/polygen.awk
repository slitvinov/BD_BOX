#! /usr/bin/awk -f
# generate an initial configuration of the polyme
BEGIN {
# polymer
    kb = 1.38e-23
    Rconst = 8.3144 # J * mol^-1 * K^-1
    jInkcal = 4184
    T = 298.15
    
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
    printf("(polygen.awk) dy: %e\n", dy) > "/dev/stderr"

    # the central charge
    Q=0.0
    # Lennard-Jones well depth
    eps= Rconst/jInkcal*T
    printf("eps: %f\n", eps) > "/dev/stderr"

    # mass
    mass=1
    # the equilibrium bond length
    r0=0.0
    # the maximum bond length (R is doubled  LJ, see documentation)
    rmax=1.5*R
    # the force constant
    H= 123 * (Rconst*T/jInkcal) / (rmax*rmax)
    printf("H: %f\n", H) > "/dev/stderr"
    # the hydrodynamic radius of the subunit
    sigma= rmax/3.0
    printf("sigma: %f\n", sigma) > "/dev/stderr"

    x0=0.0
    y0=-N/2.0*dy
    z0=0.0

    for (id=1; id<N+1; id++) {
	printf("sub DNA %i %e %e %e %e %e %e %e %e\n", 
	       id, x0, y0, z0, sigma, Q, R, eps, mass)
	y0+=dy
    }

    for (id=1; id<N; id++) {
	printf("bond %i %i %e %e %e\n",
	       id, id+1, r0, rmax, H);
    }
}
