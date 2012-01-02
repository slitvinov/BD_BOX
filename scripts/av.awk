#! /usr/bin/awk -f

# average data in several files
BEGIN {
    if (length(icol)==0) {
	icol=1
    }
}

NF{
    rest[FNR, icol] += $(icol)
    n[FNR]+=1

    if (FNR>maxfnr) {
	maxfnr=FNR
    }

    $(icol) = "" 
    time[FNR] = $0
}

END {
    for (q=1; q<maxfnr+1; ++q) {
	print time[q], rest[q, icol]/n[q]
    }
}