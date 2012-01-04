#! /usr/bin/awk -f

# Transform XYZ file to punto format
# http://sourceforge.net/projects/punto/
# http://en.wikipedia.org/wiki/XYZ_file_format

# get centro of the mass (xc, yc, zc)
function get_cm(                i) {
    xc=yc=zc=0
    for (i=1; i<ib+1; i++) {
	xc+=x[i]
	yc+=y[i]
	zc+=z[i]
    }
    xc/=ib
    yc/=ib
    zc/=ib
}

# shift the center of mass
function move_cm(               i) {
    for (i=1; i<ib+1; i++) {
	x[i]= x[i] - xc
	y[i]= y[i] - yc
	z[i]= z[i] - zc
    }
}

# get the first Rouse vector
function get_fr(                i,N,pi,fcos) {
    N=ib
    xr=yr=zr
    pi=3.141592653589793
    for (i=1; i<N+1; i++) {
	#fcos= cos (pi* (2*i+1) / (2*N) )
	fcos= cos (pi* (i-0.5) / N )
	xr+= x[i] * fcos
	yr+= y[i] * fcos
	zr+= z[i] * fcos
    }
    xr/=N
    yr/=N
    zr/=N
}

function end2endY(                         ) {
    ye = y[ib] - y[1]
}

function print_snap(            i) {
    if (e2e==1) {
	print ye
    } else if (fr==1) {
	print xr, yr, zr
    } else {
	for (i=1; i<ib+1; i++) {
	    if (i==1) {
		print x[i], y[i], z[i], col_head
	    } else if (i==ib) {
		print x[i], y[i], z[i], col_tail
	    } else {
		print x[i], y[i], z[i], col_rest
	    }
	}
	printf("\n")	
    }
} 

BEGIN {
    next_is_comment=0
    ib=0
    col_head=749
    col_tail=10
    col_rest=2
}

next_is_comment{
    # skip a commnet line
    next_is_comment=0
    next
}

NF==1{
    # <number of atoms> line
    next_is_comment=1
    if (NR>1) {
	if (cm==1) {
	    get_cm()
	    move_cm()
	}
	if (fr==1) {
	    get_fr()
	}
	if (e2e==1) {
	    end2endY()
	}
	print_snap()
    }
    ib=0
    next
}

{
    ib++
    # a line with data
    x[ib]=$2
    y[ib]=$3
    z[ib]=$4
}
