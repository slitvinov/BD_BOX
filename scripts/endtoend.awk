#! /usr/bin/awk -f

# Prints end-to-end vector
# ns: number of snapshots to skip

BEGIN {

}

function end_of_polymer() {
    print sqrt( (x-xh)^2 + (y-yh)^2 + (z-zh)^2)
    nsnap++
    npoly++
    nb=0
}

FNR==1 && NR>1 {
    # to use with several files
    printf "(endtoend.awk) newfile\n" > "/dev/stderr"
    npoly=0
}

npoly>=ns{
    
}

NF{
    nb++
    x=$1
    y=$2
    z=$3

    if (nb==1) {
	xh=x
	yh=y
	zh=z
    }
}


!NF {
    end_of_polymer()
}