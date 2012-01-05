BEGIN{
    if (length(idx)==0) {
	idx=1
    }
}

NF{
    array[N]=$(idx);
    N++
}


END  {
    for (i=0; i<N; i++) {
	sum2 += array[i]^2
    }

    for (s=0; s<N; s++) {
	sum1=0.0
	for(i=0; i<N-s; i++){
	    sum1 += (array[i+s]-array[i])^2
	} 
	print 0.5*sum1/(N-s)/(sum2/N)
    }
}
