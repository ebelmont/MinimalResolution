//unify a vector to make a given entry 1
template<typename ring>
vectors<matrix_index,ring> matrix_mem<ring>::unify(matrix_index n, vectors<matrix_index,ring> const &x){
	ring po = this->moduleOper->ringOper->inverse(this->moduleOper->component(n,x));
	return this->moduleOper->scalor_mult(po,x);
}

//ellimination a given entry
template<typename ring>
void matrix_mem<ring>::take_away(matrix_index n, vectors<matrix_index,ring> &x, vectors<matrix_index,ring> const &y){
	ring cp = this->moduleOper->ringOper->minus(this->moduleOper->component(n,x));
	auto smd = this->moduleOper->scalor_mult(cp,y);
	x = this->moduleOper->add(std::move(x),std::move(smd));
}

// row transfermations, subtract the given row to cancel the given column
template<typename ring>
void matrix_mem<ring>::row_reduction(matrix_index row, matrix_index col){
	//unify the given row
	auto uv = unify(col,data[row]);

	for(unsigned i=0; i<data.size(); ++i){
		//unify the given row
		if(i==row) data[i] = uv;
		//transform the other rows with vanhishing given column
		else take_away(col, data[i], uv);
	}
}

//Gaussian ellimination, with given rows and colums.
//
//NOTE: an earlier attempt at parallelizing this wrapped row_reduction's row
//loop in its own `#pragma omp parallel for`, spawning a fresh thread team on
//*every* pivot -- measured 55-65% slower at degree 80 (pure fork/join tax,
//see NOTES.md's 2026-08-29 entry). A second attempt (one persistent team
//for the whole call, `#pragma omp for` per pivot instead of a team spawn)
//fixed that, but on a *shared/contended* node still lost to serial: with
//only one barrier per pivot, and row_cols in the thousands per stage,
//OS scheduling latency under contention (other users' processes competing
//for the same cores) dominates, showing up as CPU-seconds far exceeding
//the serial baseline's (threads burning cycles waiting at each barrier).
//
//Fix: batch BLOCK pivots between synchronization points, same idea as
//blocked/panel LU factorization (LAPACK's dgetrf). A block's pivots are
//first fully Gauss-Jordan-reduced against EACH OTHER in a small serial
//"panel" step (touches only the ~BLOCK pivot rows, cheap) -- this leaves
//each pivot row with a zero in every OTHER pivot column of the block. That
//makes the block's combined effect on any other row order-independent: b
//sequential take_away calls using the finalized pivot rows, done back to
//back with no synchronization between them, give the exact same result as
//b separate barriered eliminations (since none of the b pivot rows touch
//each other's columns any more). So the whole block needs exactly ONE
//barrier (for the trailing `#pragma omp for`) instead of one per pivot --
//amortizing barrier/scheduling cost over BLOCK pivots at a time.
template<typename ring>
void matrix_mem<ring>::gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols){
	int nthreads = omp_get_max_threads();
	size_t total = row_cols.size();
	if(nthreads <= 1 || data.size() < 100 || total < 100){
		int i = 0;
		for(auto rc: row_cols){
			//throttled: row_cols can be in the thousands per resolution stage x
			//220 stages, and production stdout is typically redirected to a log
			//file on a network filesystem, where an unthrottled per-iteration
			//flush becomes a large number of remote-filesystem writes.
			++i;
			if(i % 500 == 0)
				std::cout << "\r" << i << "/" << total << std::flush;
			row_reduction(rc.first,rc.second);
		}
		return;
	}

	DiagTimer _dt("gaussian_parallel", nthreads);
	int n = (int)data.size();
	//block size: batches BLOCK pivots per synchronization point (see design
	//comment above). The panel step is serial and costs O(BLOCK) per pivot
	//summed over a call (BLOCK^2 work per block x total/BLOCK blocks); the
	//parallel trailing update costs O(n/nthreads) per pivot regardless of
	//BLOCK. So the panel stays cheap only while BLOCK stays well below
	//n/nthreads -- past that "crossover" a bigger BLOCK stops helping and
	//starts hurting (measured directly: at n=1572/nthreads=8, crossover
	//~196, and BLOCK=256 already lost to BLOCK=64; at n=2731, crossover
	//~341, and BLOCK=256 only tied BLOCK=64 instead of beating it -- see
	//NOTES.md's block-size tuning entry). A fixed constant can't be right
	//across stages ranging from rank ~1500 to 5000+, so scale with the
	//matrix's own size, staying a safety factor below the crossover.
	//GAUSSIAN_BLOCK env var overrides with an absolute value (tuning only).
	static const int gaussian_block_override = [](){
		const char *e = std::getenv("GAUSSIAN_BLOCK");
		return e ? std::atoi(e) : 0;
	}();
	size_t BLOCK;
	if(gaussian_block_override > 0){
		BLOCK = (size_t)gaussian_block_override;
	} else {
		int adaptive = n / (nthreads * 4);
		BLOCK = (size_t)(adaptive < 16 ? 16 : (adaptive > 512 ? 512 : adaptive));
	}
	//breakdown instrumentation (DIAG_TIMING-gated, no effect on production
	//runs): separates the call's total wall/cpu (already reported by the
	//DiagTimer above) into (1) time in the serial panel step, (2) genuine
	//per-thread CPU spent in the parallel trailing step (RUSAGE_THREAD, so
	//it can't be inflated by other threads' idle/waiting time), leaving
	//whatever's unaccounted-for as synchronization/scheduling overhead.
	bool diag_on = Diag::enabled();
	double panel_wall_total = 0.0, panel_cpu_total = 0.0;
	std::vector<double> trailing_cpu_per_thread(diag_on ? nthreads : 0, 0.0);
	std::vector<double> trailing_wall_per_thread(diag_on ? nthreads : 0, 0.0);
	//element-visit counters (not timing -- exact counts, negligible overhead
	//vs. an extra getrusage call per row): distinguishes genuinely-necessary
	//accumulate work (cp nonzero, a real contribution) from add2Dense's
	//unconditional per-pivot-row-element visits when cp is zero (no math
	//needed, but the current code has no early-out and pays the visit
	//anyway), and separately tallies the populate/de-densify merge costs
	//that have no equivalent in the old take_away-per-pivot design at all.
	std::vector<long long> accumElemsUseful(diag_on ? nthreads : 0, 0);
	std::vector<long long> accumElemsWasted(diag_on ? nthreads : 0, 0);
	std::vector<long long> accumCallsUseful(diag_on ? nthreads : 0, 0);
	std::vector<long long> accumCallsWasted(diag_on ? nthreads : 0, 0);
	std::vector<long long> populateElems(diag_on ? nthreads : 0, 0);
	std::vector<long long> dedenseElems(diag_on ? nthreads : 0, 0);
	//direct ground-truth check for whether the omp-for worksharing is
	//actually partitioning the n-row trailing loop across threads (each
	//thread gets a disjoint subset) vs. some bug making every thread
	//redundantly execute the whole loop -- counts rows *actually entered*
	//per thread, no assumption about how work should be divided. Also
	//records each thread's own reported team size once, to rule out a
	//degenerate 1-thread team silently satisfying `nthreads` elsewhere.
	std::vector<long long> rowVisits(diag_on ? nthreads : 0, 0);
	std::vector<int> observedTeamSize(diag_on ? nthreads : 0, -1);

	//per-thread dense scratch buffers for the trailing step (see that loop's
	//comment below for why): reused by every row and every block this
	//thread processes for the rest of the call. Cheap to size upfront --
	//`ring`'s default ctor (a `vectors<>`) allocates nothing, so this is
	//one allocation of bufSize*sizeof(ring) per thread, not per-element
	//allocations.
	//
	//Sized to the *true* max column index actually used, NOT to `n`
	//(=data.size(), the row count): column indices are a separate index
	//space from row indices for this matrix (confirmed the hard way -- an
	//explicit bounds check below caught column indices exceeding `n` on
	//the very first correctness test; `indj`'s columns come from an
	//embedding map's own index space, which is not guaranteed to match its
	//row count). Scanning every row plus every pivot column once, up
	//front, is O(n + pivots) -- negligible next to the O(n*pivots) main
	//work below.
	matrix_index maxCol = 0;
	for(auto &row : data)
		if(!row.dataArray.empty())
			maxCol = std::max(maxCol, row.dataArray.back().ind);
	for(auto &rc : row_cols)
		maxCol = std::max(maxCol, rc.second);
	size_t bufSize = (size_t)maxCol + 1;
	std::vector<std::vector<ring>> denseBuf(nthreads);
	for(int t=0; t<nthreads; ++t) denseBuf[t].resize(bufSize);
	//per-thread scratch for a single row's own (already-sorted) column
	//indices, merged on the fly against this row's own touched-column set
	//(see touchedBuf below, and deDenseFromSortedUnion) when converting
	//back to sparse -- `bufSize` (the column space) can be several times
	//larger than any single row's density, so that conversion must only
	//visit indices that could actually be nonzero, not scan the whole
	//buffer.
	std::vector<std::vector<matrix_index>> rowIdxBuf(nthreads);
	//per-thread reused scratch for the *per-row* touched-column set (see
	//add2DenseTracked in modules/2.h): collected while applying this
	//block's pivots to one row, sized to that row's own actual footprint
	//(sum of the densities of only the pivots with a nonzero coefficient
	//against it) rather than the whole block's shared union. A prior
	//version of this code used one blockTouched union (computed once per
	//block, shared by every row) here instead -- correct, but every row
	//paid O(blockTouched.size()) in deDenseFromSortedUnion regardless of
	//how much of that union actually applied to it; this per-row
	//alternative was measured to fix exactly that. Reused (cleared, not
	//reallocated) across rows and blocks by each thread; grows by
	//amortized doubling on first use and then keeps its capacity.
	std::vector<std::vector<matrix_index>> touchedBuf(nthreads);

	#pragma omp parallel num_threads(nthreads)
	{
		int tid = omp_get_thread_num();
		if(diag_on) observedTeamSize[tid] = omp_get_num_threads();
		std::vector<ring> &buf = denseBuf[tid];
		std::vector<matrix_index> &rowIdx = rowIdxBuf[tid];
		std::vector<matrix_index> &touched = touchedBuf[tid];
		for(size_t base=0; base<total; base+=BLOCK){
			size_t blockEnd = std::min(base+BLOCK, total);
			//panel step: serial, touches only this block's own pivot rows --
			//identical math to row_reduction but scoped to the block instead
			//of the whole matrix (see comment above for why that's sufficient)
			#pragma omp single
			{
				//one print per block (<=BLOCK pivots), not per pivot -- already
				//throttled relative to the original per-500-pivot cadence
				std::cout << "\r" << blockEnd << "/" << total << std::flush;
				std::chrono::steady_clock::time_point pt0;
				double pc0 = 0.0;
				if(diag_on){ pt0 = std::chrono::steady_clock::now(); pc0 = Diag::thread_cpu_seconds(); }
				for(size_t k=base; k<blockEnd; ++k){
					matrix_index row_k = row_cols[k].first, col_k = row_cols[k].second;
					data[row_k] = unify(col_k, data[row_k]);
					for(size_t j=base; j<blockEnd; ++j){
						if(j==k) continue;
						take_away(col_k, data[row_cols[j].first], data[row_k]);
					}
				}
				if(diag_on){
					panel_wall_total += std::chrono::duration<double>(std::chrono::steady_clock::now()-pt0).count();
					panel_cpu_total += Diag::thread_cpu_seconds()-pc0;
				}
			}
			//trailing update: parallel over every row NOT among this block's
			//own pivots (includes earlier/later blocks' pivot rows too, same
			//as row_reduction's original "all i != row" semantics).
			//schedule(runtime) (tunable via OMP_SCHEDULE, e.g.
			//"OMP_SCHEDULE=static,1" for a cyclic/round-robin assignment
			//instead of the default's one-contiguous-chunk-per-thread) --
			//plain schedule(static) hands each thread a contiguous run of
			//row indices, so if row "difficulty" (nonzero count) has any
			//spatial correlation with index (plausible: rows are ordered by
			//generator/degree, and nearby generators can share term-count
			//profiles), one thread can get handed a run of all-dense rows
			//while another gets all-sparse ones. A cyclic chunk of 1 hedges
			//against that by spreading any such correlated run across every
			//thread instead of concentrating it in one.
			double tc0 = diag_on ? Diag::thread_cpu_seconds() : 0.0;
			std::chrono::steady_clock::time_point tw0;
			if(diag_on) tw0 = std::chrono::steady_clock::now();
			//NOTE on the dense-buffer accumulation below: profiling (see
			//NOTES.md's 2026-08-30 entry) found the dominant cost of this
			//loop was NOT synchronization or arithmetic, but `take_away`
			//building a fresh, freshly-allocated full-row result vector on
			//EVERY pivot application (once per (row,pivot) pair) -- most of
			//that cost was in trivially-O(1) vector/term construction and
			//move machinery, made expensive by 8 threads concurrently
			//first-touching cold memory at that volume. Fix: accumulate
			//this block's BLOCK pivots into one reused dense buffer per
			//thread (populate once from data[r], `add2Dense` per pivot with
			//zero allocation, convert back to sparse once) instead of
			//allocating a fresh sparse result vector per pivot.
			#pragma omp for schedule(runtime)
			for(int r=0; r<n; ++r){
				bool isBlockPivot = false;
				for(size_t k=base; k<blockEnd && !isBlockPivot; ++k)
					if(row_cols[k].first == (matrix_index)r) isBlockPivot = true;
				if(isBlockPivot) continue;
				//snapshot this row's own column indices (already sorted,
				//since dataArray always is) before data[r] gets overwritten
				//below -- deDenseFromSortedUnion merges this on the fly
				//against `touched` (this row's own actually-touched column
				//set, collected below as this block's pivots are applied --
				//see touchedBuf's declaration comment for why per-row
				//instead of per-block matters for cost, not correctness).
				rowIdx.clear();
				touched.clear();
				for(auto &t : data[r].dataArray){
					//defensive: bufSize is computed from the same data up
					//front, so this should never fire -- checked explicitly
					//since a violation would otherwise be a silent OOB
					//write, not just a wrong answer.
					if(t.ind >= (matrix_index)bufSize){
						std::cerr << "FATAL: dense buffer column index " << t.ind
						          << " >= bufSize=" << bufSize << " in gaussian() trailing step" << std::endl;
						std::abort();
					}
					buf[t.ind] = std::move(t.coeficient);
					rowIdx.push_back(t.ind);
				}
				if(diag_on){ populateElems[tid] += (long long)rowIdx.size(); rowVisits[tid] += 1; }
				for(size_t k=base; k<blockEnd; ++k){
					matrix_index row_k = row_cols[k].first, col_k = row_cols[k].second;
					if(col_k >= (matrix_index)bufSize){
						std::cerr << "FATAL: pivot column " << col_k << " >= bufSize=" << bufSize
						          << " in gaussian() trailing step" << std::endl;
						std::abort();
					}
					ring cp = this->moduleOper->ringOper->minus(buf[col_k]);
					if(diag_on){
						long long pivotDensity = (long long)data[row_k].dataArray.size();
						if(this->moduleOper->ringOper->isZero(cp)){
							accumElemsWasted[tid] += pivotDensity;
							accumCallsWasted[tid] += 1;
						} else {
							accumElemsUseful[tid] += pivotDensity;
							accumCallsUseful[tid] += 1;
						}
					}
					data[row_k].add2DenseTracked(cp, buf, this->moduleOper->ringOper, touched);
				}
				//sort+dedup this row's actually-touched columns (typically
				//much smaller than the old whole-block union) before
				//merging with rowIdx -- deDenseFromSortedUnion requires
				//both lists sorted, and tolerates but doesn't need
				//dedup'd input; dedup'ing here keeps its merge loop from
				//doing extra no-op work on repeated indices (the same
				//column can be touched by more than one of this block's
				//pivots).
				std::sort(touched.begin(), touched.end());
				touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
				if(diag_on) dedenseElems[tid] += (long long)(rowIdx.size() + touched.size());
				data[r].deDenseFromSortedUnion(buf, this->moduleOper->ringOper, rowIdx, touched);
			}
			if(diag_on){
				//taken right after the `for`'s own implicit barrier -- so this
				//thread's wall delta includes time spent waiting for slower
				//threads to finish their share, distinguishing genuine load
				//imbalance (large spread across threads) from pure
				//synchronization/wake latency (all threads similar, but all
				//longer than their own cpu time would predict)
				trailing_wall_per_thread[tid] += std::chrono::duration<double>(std::chrono::steady_clock::now()-tw0).count();
				trailing_cpu_per_thread[tid] += Diag::thread_cpu_seconds()-tc0;
			}
		}
	}
	if(diag_on){
		double trailing_cpu_total = 0.0;
		double trailing_wall_min = trailing_wall_per_thread.empty() ? 0.0 : trailing_wall_per_thread[0];
		double trailing_wall_max = trailing_wall_min;
		for(double v : trailing_cpu_per_thread) trailing_cpu_total += v;
		for(double v : trailing_wall_per_thread){
			if(v < trailing_wall_min) trailing_wall_min = v;
			if(v > trailing_wall_max) trailing_wall_max = v;
		}
		//trailing_wall_max is what a thread that's never waiting would need --
		//i.e. the load-imbalance floor. If it's close to trailing_cpu_total/
		//nthreads, work is balanced and the remaining gap to the top-level
		//DiagTimer's wall is genuine sync/wake overhead; if trailing_wall_max
		//is itself much larger than trailing_wall_min, that's load imbalance,
		//not synchronization cost, and schedule(static) is the thing to fix
		std::cerr << "DIAG gaussian_breakdown panel_wall=" << panel_wall_total
		          << " panel_cpu=" << panel_cpu_total
		          << " trailing_cpu_sum=" << trailing_cpu_total
		          << " trailing_wall_est=" << (trailing_cpu_total/nthreads)
		          << " trailing_wall_min=" << trailing_wall_min
		          << " trailing_wall_max=" << trailing_wall_max
		          << " nthreads=" << nthreads
		          << std::endl;
		long long su=0, sw=0, cu=0, cw=0, pop=0, ded=0, rv=0;
		for(int t=0;t<nthreads;++t){
			su+=accumElemsUseful[t]; sw+=accumElemsWasted[t];
			cu+=accumCallsUseful[t]; cw+=accumCallsWasted[t];
			pop+=populateElems[t]; ded+=dedenseElems[t];
			rv+=rowVisits[t];
		}
		std::cerr << "DIAG gaussian_workvisits accum_elems_useful=" << su
		          << " accum_elems_wasted=" << sw
		          << " accum_calls_useful=" << cu
		          << " accum_calls_wasted=" << cw
		          << " populate_elems=" << pop
		          << " dedense_elems=" << ded
		          << " row_visits_total=" << rv
		          << " expected_row_visits_if_partitioned=" << ((total/BLOCK + 1) * (long long)n)
		          << std::endl;
		std::cerr << "DIAG gaussian_teamsizes";
		for(int t=0;t<nthreads;++t) std::cerr << " t" << t << "=" << observedTeamSize[t];
		std::cerr << std::endl;
	}
}
