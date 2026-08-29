#include"matrices_mem.h"
#include<omp.h>
#include<algorithm>

//compute the quotient comodule. Parallelized across the surviving
//generators (inverse_ind): the per-generator work (quot->maps_to(...))
//is independent, but two hazards rule out a naive parallel-for --
//
// 1) `quot` may be a matrix_file (stream-backed; see matrices_stream.h),
//    whose find()/maps_to() go through one con_streams object with a
//    single shared seek position. con_streams::read() does hold its own
//    mutex around seek+read, so concurrent calls aren't a *data race*,
//    but they'd fully serialize on that lock -- worse than useless once
//    parallelized, since it adds lock overhead on top of the same
//    effective serial I/O. Fix: materialize `quot` into a matrix_mem
//    once, up front (same pattern as quotient_p below), so the parallel
//    region's reads are plain in-memory lookups.
// 2) `result.coaction_matrix` (the write target) can likewise be a
//    matrix_file (e.g. BP_init.h's resolution comodule) -- writes must
//    stay serial, not be raced from multiple threads.
//
// X->coaction(i) was checked too (traced to cofree_comodule::coaction,
// hopf_algebroid/13.h): it's a const function of i, but for the real
// BP_*BP algebroid it bottoms out in BP_Op::delta(n) = delta_table->find(n),
// and delta_table is itself a matrix_file (BP_init.h's `delta_matrix`) --
// so it gets the same stream hazard as (1) and is gathered serially below
// rather than raced inside the parallel region.
//
// (Z3_Op::power_p's cache, unsafe to grow concurrently -- see Z3.cpp and
// the "Parallelize the bucketed Curtis table search" commit -- was also
// checked here: the only power_p/divide_power_p call site in BP.cpp is
// line ~247, inside theta-class construction, not reachable from
// delta()/etaR()/right_scalor_mult()/multiply() on this path. No pre-warm
// needed.)
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::quotient(const CoModule<algebroid,degree_type> *X, matrix<ring> *quot, std::vector<matrix_index> inverse_ind, comodule_generic<algebroid,degree_type> &result){
	//get the rank of the quotient comodule
	result.base_module.rank = inverse_ind.size();

	//compute the degrees of the generators of the quotient
	result.base_module.degree.resize(inverse_ind.size());
	for(unsigned i=0; i<inverse_ind.size(); ++i)
		result.base_module.degree[i] = X->degree(inverse_ind[i]);
	std::function<vectors<matrix_index,algebroid>(const algebroid&, const vectors<matrix_index,ring>&)> right_mult = [this] (const algebroid& A, const vectors<matrix_index,ring> &V){
		return right_scalor_mult(A,V); };

	unsigned n = (unsigned)inverse_ind.size();
	AbGroupOp<vectors<matrix_index,algebroid>> *adop = algebroidModuleOper;

	//materialize quot (hazard 1 above)
	matrix_mem<ring> quot_mem;
	quot_mem.construct(quot);

	//serially gather each generator's coaction (hazard noted above); this
	//also hands us term counts for LPT load-balancing below for free
	std::cout << "computing the coactions...\n" << std::flush;
	std::vector<vectors<matrix_index,algebroid>> coas(n);
	std::vector<size_t> term_count(n);
	for(unsigned i=0; i<n; ++i){
		coas[i] = X->coaction(inverse_ind[i]);
		term_count[i] = coas[i].size();
	}

	//LPT partition: sort indices descending by coaction term count (a
	//direct proxy for maps_to's per-generator cost), greedily assign to
	//the currently lightest thread. This sort+assign runs once per
	//resolvor_modeled call over `n` surviving generators, not once per
	//coaction term -- see NOTES.md's "quotient()'s LPT load-balancing"
	//entry for why that's not the earlier row-level-OpenMP mistake.
	int nthreads = std::max(1, std::min(omp_get_max_threads(), std::max((int)n,1)));
	std::vector<int> order(n);
	for(unsigned i=0; i<n; ++i) order[i]=i;
	std::sort(order.begin(), order.end(), [&](int a, int b){
		return term_count[a] > term_count[b]; });
	std::vector<std::vector<int>> assign(nthreads);
	std::vector<size_t> load(nthreads, 0);
	for(int idx : order){
		int t = (int)(std::min_element(load.begin(), load.end()) - load.begin());
		assign[t].push_back(idx);
		load[t] += term_count[idx];
	}

	//parallel compute: each thread writes only the positions assigned to
	//it in `results`, indexed by position within inverse_ind -- no shared
	//mutation races
	std::vector<vectors<matrix_index,algebroid>> results(n);
	#pragma omp parallel num_threads(nthreads)
	{
		int tid = omp_get_thread_num();
		for(int i : assign[tid]){
			//throttled -- see matrices_mem/4.h::gaussian
			if(i % 500 == 0)
				std::cout << "\r" << i << "/" << n << std::flush;
			results[i] = quot_mem.maps_to(coas[i], right_mult, adop);
		}
	}

	//serial writeback (hazard 2 above)
	std::function<vectors<matrix_index,algebroid>(int)> rows = [&results](int i){
		return results[i]; };
	result.coaction_matrix->construct((int)n, rows);
}

//compute the quotient comodule
template<typename ring, typename algebroid>
template<typename degree_type>
void Hopf_Algebroid<ring,algebroid>::quotient_p(const CoModule<algebroid,degree_type> *X, matrix<ring> *quot0, std::vector<matrix_index> inverse_ind, comodule_generic<algebroid,degree_type> &result){
	//load the quot matrix
	matrix_mem<ring> quot1;
	matrix<ring> *quot = &quot1;
	quot->construct(quot0);
	
	//get the rank of the quotient comodule
	result.base_module.rank = inverse_ind.size();
	
	//compute the degrees of the generators of the quotient
	result.base_module.degree.resize(inverse_ind.size());
	for(unsigned i=0; i<inverse_ind.size(); ++i)
		result.base_module.degree[i] = X->degree(inverse_ind[i]);
	std::function<vectors<matrix_index,algebroid>(const algebroid&, const vectors<matrix_index,ring>&)> right_mult = [this] (const algebroid& A, const vectors<matrix_index,ring> &V){
		return right_scalor_mult(A,V); };
		
		//compose the quotien matrix and the coaction matrix. Note that we need to first apply the right unit to the quotient matrix
		std::function<vectors<matrix_index,algebroid>(int)> coactor = [X,quot,&inverse_ind,right_mult,this] (int i){
			vectors<matrix_index,algebroid> coa = X->coaction(inverse_ind[i]);
			AbGroupOp<vectors<matrix_index,algebroid>> *adop = algebroidModuleOper;
			//throttled -- see matrices_mem/4.h::gaussian
			if(i % 500 == 0)
				std::cout << "\r" << i << "/" << inverse_ind.size() << std::flush;
			return quot->maps_to_p(coa,right_mult, adop);
		};
		std::cout << "computing the coactions...\n" << std::flush;
		result.coaction_matrix->construct_parallel(inverse_ind.size(), coactor);
}

//scalor multiplication using the right unit
template<typename ring, typename algebroid>
vectors<matrix_index, algebroid> Hopf_Algebroid<ring,algebroid>::right_scalor_mult (algebroid const &A, vectors<matrix_index, ring> const &V) {
	if(algebroidRingOper->isZero(A))
		return algebroidModuleOper->zero();
	std::function<algebroid(ring const&)> rule = [this,&A] (ring const &r) { return algebroidRingOper->multiply(A,etaR(r)); };
	return moduleOper->termwise_operation(rule, V);
}
