# get msd of the center of mass

function r2(x, y, z) {
    return x^2+y^2+z^2
}

NF{
    x[N]=$1
    y[N]=$2
    z[N]=$3
    N++
}


END  {
    for (s=0; s<N; s++) {
	sum1=0.0
	for(i=0; i<N-s; i++){
	    sum1 += r2(x[i+s]-x[i], y[i+s]-y[i], z[i+s]-z[i])
	} 
	print 0.5*sum1/(N-s)
    }
}
