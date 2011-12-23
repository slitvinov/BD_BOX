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

function move_cm(               i) {
    for (i=1; i<ib+1; i++) {
	x[i]= x[i] - xc
	y[i]= y[i] - yc
	z[i]= z[i] - zc
    }
}

function print_snap(            i) {
    if (cm==1) {
	get_cm()
	move_cm()
    }
    for (i=1; i<ib+1; i++) {
	if (i==1) {
	    print x[i], y[i], z[i], col_head
	} else if (i==ib) {
	    print x[i], y[i], z[i], col_tail
	} else {
	    print x[i], y[i], z[i], col_rest
	}
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
	print_snap()
	printf("\n")
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
