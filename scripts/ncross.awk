#! /usr/bin/awk -f

# number of zero line crossing

function sgn(x) {
    if (x>0) {
	return(1)
    } else {
	return(-1)
    }
}

NR==1 {
    # previous sign
    psgn= sgn($1)
    next
}

NF{ 
    csgn = sgn($1)
    if (csgn != psgn) {
	# this is crossing
	nc++
    }
    psgn=csgn
    n++
}

END {
    # average time between crossing
    print n/nc
}


