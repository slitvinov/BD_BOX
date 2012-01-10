#! /usr/bin/awk -f

# Get average polymer extension in X, Y, and Z directions
# ns: number of snapshots to skip

BEGIN {
   maxx=maxy=maxz=-1e19
   minx=miny=minz=+1e19
}

function end_of_polymer() {
    av_extx+= (maxx - minx)^2
    av_exty+= (maxy - miny)^2
    av_extz+= (maxz - minz)^2

    if (debug) {
	print maxx, maxy, maxz, minx, miny, minz
    }
    maxx=maxy=maxz=-1e19
    minx=miny=minz=+1e19

    nsnap++
    npoly++
}

FNR==1 && NR>1 {
    # to use with several files
    printf "(ext.awk) newfile\n" > "/dev/stderr"
    npoly=0
}

npoly>=ns{
    
}

NF{
    x=$1
    y=$2
    z=$3

    if (x>maxx) maxx=x
    if (x<minx) minx=x

    if (y>maxy) maxy=y
    if (y<miny) miny=y

    if (z>maxz) maxz=z
    if (z<minz) minz=z

}


!NF {
    end_of_polymer()
}

END {
    print av_extx/nsnap,
	av_exty/nsnap,
	av_extz/nsnap
}