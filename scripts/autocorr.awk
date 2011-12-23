NF{ 
    array[N]=$1; 
    N++
}


END  {
    for(i=0; i<N; i++){
	sum += array[i];
    } 
    mean = sum/N;
    for (i=0; i<N; i++) {
	sum2 += (array[i]-mean)^2;
    }

    for (s=0; s<N; s++) {
	sum1=0.0
	for(i=0; i<N-s; i++){
	    sum1 += ((array[i+s]-mean)*(array[i]-mean));
	} 
	print sum1/sum2;
    }
	
}