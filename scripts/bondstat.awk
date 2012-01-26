#! /usr/bin/awk -f

# get bond statistics from
# XYZ file

function fabs(x) {
    if (x>0) {
	return x 
    } else {
	return -x
    }
}


function print_snap(            i) {
    for (i=2; i<ib+1; i++) {
	print sqrt((x[i]-x[i-1])^2 + (y[i]-y[i-1])^2 + (z[i]-z[i-1])^2)
    }
} 

BEGIN {
    next_is_comment=0
    ib=0
    # colors for punto
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
