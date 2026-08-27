//compute the sum of a array of elements
template<typename return_type>
return_type sum(std::function<return_type(int)> summands, int start, int end, std::function<return_type(return_type&&, return_type&&)> add){
	if(end-start==1) return summands(start);
	return_type s1,s2;
	int mid = (start+end)/2;
	//avoid ever entering the OpenMP runtime on the serial path -- an omp
	//"if()" clause still dispatches through GOMP_parallel even when false,
	//and sum()'s recursion hits this every level, so that dispatch overhead
	//dominates real work (measured: 262% CPU of 24 requested, System time
	//~4x User time, 10M+ context switches for a 130s run at degree 80)
	if(omp_site_enabled(0) && !omp_in_parallel()){
		#pragma omp parallel sections
		{
			#pragma omp section
			{ s1 = sum(summands, start,mid, add); }
			#pragma omp section
			{ s2 = sum(summands, mid,end, add); }
		}
	}
	else {
		s1 = sum(summands, start,mid, add);
		s2 = sum(summands, mid,end, add);
	}
	return add(std::move(s1), std::move(s2));
}

//the sum of an array of elements in an abelian group
template<typename A>
A sum(std::function<A(int)> summands, int start, int end, AbGroupOp<A> *abOper){
	std::function<A(A&&, A&&)> adding = [abOper] (A &&x, A &&y){ 
		return abOper->add(std::move(x),std::move(y)); };
	if(end<=start) return abOper->zero();
	return sum(summands, start, end, adding);
}

//the sum of an array of vectors
template<typename index, typename R>
vectors<index,R>  ModuleOp<index,R>::sum(std::function<vectors<index,R>(int)> summands, int start, int end){
	return ::sum(summands, start, end, this);
}
