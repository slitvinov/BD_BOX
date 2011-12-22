#! /usr/bin/awk -f

# Transform XYZ file to punto format
# http://sourceforge.net/projects/punto/
# http://en.wikipedia.org/wiki/XYZ_file_format

BEGIN {
    next_is_comment=0
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
	printf("\n")
    }
    next
}

{
    # a line with data
    # we remove tag and printf the reset
    x=$2
    y=$3
    z=$4
    print x, y, z
}
