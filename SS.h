//SS.h
#pragma once
#include"matrices.h"
#include<memory>
#include<algorithm>
#include<omp.h>

//class for the entries in a spectral sequence
template<typename cycle_name, typename ring>
class SS_entry{
public:
	//type of the full cycles
	typedef vectors<matrix_index, ring> value_type;
	
	//leading term of an entry
	cycle_name tag, cycle;
	//full cycle and full tags
	value_type full_tag, full_cycle;
};
    
//the class for a spectral sequence
template<typename cycle_name, typename ring>
class SS_table : public std::vector<SS_entry<cycle_name, ring>>{
public:
	//operations on the vectors
	ModuleOp<matrix_index, ring> *Modop;
        
	//the index of tags and cycles
	std::map<cycle_name, unsigned> tag_index, cycle_index;
	
	//the filtration function
	virtual int filtration(cycle_name)=0;

	//the (internal, t-) degree of a cycle name -- used together with
	//bucket_rows() below to partition one homological degree's search by
	//degree (the differential is degree-preserving, so degree buckets are
	//provably independent -- see algNov.cpp's bucket_by_degree)
	virtual int degree(cycle_name)=0;
	
	//the length of the differential
	int diff_length(SS_entry<cycle_name, ring> et){
		return filtration(et.cycle) - filtration(et.tag);
	}
	
	//transform a term into a cycle name
	virtual cycle_name naming(const ring &coeficient, matrix_index ind)=0;
	
	//computing the leading term of a cycle vector
	virtual unsigned leading_term(typename SS_entry<cycle_name, ring>::value_type&)=0;

	//inset an entry
	void insert(SS_entry<cycle_name,ring> et){
		auto i = this->size();
		//add the new entry in the end
		this->push_back(et);
		//update the indeces
		tag_index.emplace(et.tag,i);
		cycle_index.emplace(et.cycle,i);
	}
        
	//types of cycles
	enum cycle_type {Boundary, Untagged, NonFound};
        
	//the default invalid entry
	int const Invalid;
	virtual cycle_name invalid()=0;
	
	//the constructors
	SS_table(ModuleOp<matrix_index, ring> *modoper) : Invalid(-65537){
		Modop = modoper;	}

	//virtual: make_table's bucketed temp tables are owned via a base-class
	//(SS_table<cycle_name,ring>*) unique_ptr but are actually algNov_table/
	//Boc_table instances reached through virtual inheritance -- deleting
	//through a non-virtual base destructor there is undefined behavior
	//(observed in practice as a heap-corrupting "free(): invalid pointer"
	//crash), not merely a leak, since the delete needs the vtable to find
	//the true most-derived object address.
	virtual ~SS_table(){}
		
	//check the entry is a tag or not
	virtual bool tagged(const SS_entry<cycle_name,ring>&)=0;
	//check if an entry is valid
	virtual bool valid(const cycle_name&)=0;
	
	//the scalar s such that adding s*matched to a vector cancels its leading
	//term exactly -- i.e. s*component(ind,matched) == -coef. Unscaled
	//subtraction (just adding -matched) only cancels correctly when the two
	//leading coefficients are exactly equal, true by accident at p=2 (F2 has
	//a single nonzero residue) but not at odd primes. The two coefficients
	//share the same p-adic valuation (guaranteed by the caller: they matched
	//via cycle_name, which encodes valuation) but are not necessarily units
	//themselves (e.g. a v0-power tag has coefficient p^j, valuation j) --
	//divide out the common p^valuation factor first so inverse() sees a
	//genuine unit. NOTE: matrices/9.h and matrices_mem/4.h do NOT actually
	//use this pattern yet (they call inverse() raw) despite an earlier
	//version of this comment claiming otherwise -- see NOTES.md.
	ring cancel_scalar(matrix_index ind, const ring &coef, typename SS_entry<cycle_name, ring>::value_type &matched){
		auto matched_coef = Modop->component(ind, matched);
		unsigned val = Modop->ringOper->valuation(coef);
		auto c_unit = Modop->ringOper->divide(coef, val);
		auto m_unit = Modop->ringOper->divide(matched_coef, val);
		//the shared-valuation guarantee above is only a documented invariant,
		//not an enforced one. If it's ever violated, m_unit still carries a
		//factor of p and inverse() below is about to throw -- report exactly
		//what broke first. This check and the print are only ever reached on
		//that failure path: invertible() is the same test inverse() already
		//runs internally, so it's essentially free on the billions-of-calls
		//success path, and valuation(matched_coef) is only recomputed here,
		//not on every call.
		if(!Modop->ringOper->invertible(m_unit)){
			std::cerr << "cancel_scalar: valuation mismatch -- ind=" << ind
				<< " coef=" << Modop->ringOper->output(coef)
				<< " (val=" << val << ")"
				<< " matched_coef=" << Modop->ringOper->output(matched_coef)
				<< " (val=" << Modop->ringOper->valuation(matched_coef) << ")"
				<< "\n" << std::flush;
		}
		return Modop->ringOper->minus(Modop->ringOper->multiply(c_unit, Modop->ringOper->inverse(m_unit)));
	}

	//The residue variant of cancel_scalar, used by name_of_cycle().
	//
	//naming() records only a coefficient's VALUATION (as a v0-power) and
	//discards its unit, so the named basis is the associated graded of the
	//algebraic Novikov filtration and a named term's coefficient lives in
	//F_p, not in the full ring. Cancelling the whole leading coefficient
	//(cancel_scalar) is correct linear algebra over Z/p^k but wrong for
	//that basis: it removes the unit's entire value instead of its residue
	//mod p, so the F_p coefficient is lost and the remainder -- which
	//belongs to higher v0-powers and should resurface as its own named
	//term -- gets cancelled away with it.
	//
	//Subtracting r*matched, where r is the ratio's residue mod p, instead
	//leaves a strictly higher valuation, exactly as repeated unscaled
	//subtraction used to do implicitly before cancel_scalar was introduced
	//(the multiplicity that encoded the coefficient in the output files was
	//an emergent property of that loop, never a documented one).
	//
	//Returns {scalar to add, r}: the caller adds scalar*matched and records
	//the leading name r times, r being the term's F_p coefficient.
	std::pair<ring, unsigned> cancel_residue(matrix_index ind, const ring &coef, typename SS_entry<cycle_name, ring>::value_type &matched){
		auto matched_coef = Modop->component(ind, matched);
		unsigned val = Modop->ringOper->valuation(coef);
		auto c_unit = Modop->ringOper->divide(coef, val);
		auto m_unit = Modop->ringOper->divide(matched_coef, val);
		if(!Modop->ringOper->invertible(m_unit)){
			std::cerr << "cancel_residue: valuation mismatch -- ind=" << ind
				<< " coef=" << Modop->ringOper->output(coef)
				<< " (val=" << val << ")"
				<< " matched_coef=" << Modop->ringOper->output(matched_coef)
				<< " (val=" << Modop->ringOper->valuation(matched_coef) << ")"
				<< "\n" << std::flush;
		}
		auto ratio = Modop->ringOper->multiply(c_unit, Modop->ringOper->inverse(m_unit));
		//find r with valuation(ratio - r) >= 1, i.e. r == ratio mod p. The
		//bound is a guard, not a real limit: r < p always.
		for(unsigned r=1; r<64; ++r){
			auto probe = Modop->ringOper->unit((int)r);
			auto diff = Modop->ringOper->add(ratio, Modop->ringOper->minus(probe));
			if(Modop->ringOper->isZero(diff) || Modop->ringOper->valuation(diff) >= 1)
				return std::make_pair(Modop->ringOper->minus(probe), r);
		}
		std::cerr << "cancel_residue: no residue r<64 with valuation(ratio-r)>=1 "
			"-- ratio=" << Modop->ringOper->output(ratio) << " is not a unit of a "
			"small-residue-field ring, which this naming scheme assumes\n" << std::flush;
		std::abort();
	}

	//simplify a cycle by subtraction boundaries using the data in the table.
	//Iterative, not recursive: this used to call itself once per leading
	//term canceled, but several locals here (led_name, it, scalar, nc, nt)
	//have non-trivial destructors that must run after the "recursive" call
	//returns, so the compiler can't tail-call-optimize it even at -O2 --
	//each step was a real stack frame. A long enough simplification chain
	//(observed 8000+ steps at degree 380) overflows the stack and
	//segfaults partway through, at whatever function happens to touch the
	//now-unmapped stack memory next (seen both inside table construction
	//and, confusingly, inside malloc()) -- not a data race, despite the
	//parallel code nearby. See NOTES.md.
	cycle_type simplify(typename SS_entry<cycle_name, ring>::value_type& cyc, typename SS_entry<cycle_name, ring>::value_type& htpy, int pric, int depth=0){
		while(true){
			//if cycle is zero, then it is a boundary
			if(Modop->isZero(cyc)) return Boundary;

			//compute the leading term
			unsigned led = leading_term(cyc);
			cycle_name led_name = naming(cyc.dataArray[led].coeficient, cyc.dataArray[led].ind);

			//if the filtration exceeds the required pricision, then regard it as a boundary
			if(filtration(led_name)>pric)
				return Boundary;

			//if not trivial, try to find it in the table
			auto it = cycle_index.find(led_name);
			if(it==cycle_index.end()) return NonFound;

			//if the searched entry is not tagged, then it is a nontrivial cycle
			auto pas = it->second;
			if(this->at(pas).tag.second==Invalid) return Untagged;

			//then we subtract the tagged entry from the cycle, scaled to exactly
			//cancel the leading term (see cancel_scalar)
			auto scalar = cancel_scalar(cyc.dataArray[led].ind, cyc.dataArray[led].coeficient, this->at(pas).full_cycle);
			auto nc = Modop->scalor_mult(scalar, this->at(pas).full_cycle);
			auto nt = Modop->scalor_mult(scalar, this->at(pas).full_tag);

			cyc = Modop->add(std::move(cyc),std::move(nc));
			htpy = Modop->add(std::move(htpy),std::move(nt));

			//then simplify the subtracted cycle
			++depth;
		}
	}
        
	//clear the table
	void clearing(){
		tag_index.clear();
		cycle_index.clear();
		this->clear();
	}
       
	//output entry
	virtual std::pair<std::pair<int,int>,string> output(SS_entry<cycle_name, ring>&, int)=0;
	
	//output of cycle names
	virtual string output(cycle_name,int)=0;
	
	//IO of entries
	virtual void save(SS_entry<cycle_name, ring>&, std::iostream&)=0;
	virtual void load(SS_entry<cycle_name, ring>&, std::iostream&)=0;
	
	//IO operations
	void save(std::iostream& fl){
		unsigned ls = this->size();
		fl.write((char*)&ls,sizeof(unsigned));
		for(auto &tm: *this){
			save(tm,fl);
		}
	}
	
	//IO operations
	void load(std::iostream& fl){
		clearing();
		unsigned ls;
		fl.read((char*)&ls, sizeof(unsigned));
		SS_entry<cycle_name, ring> tm;
		for(unsigned i=0; i<ls; ++i){
			load(tm,fl);
			insert(tm);
		}
	}
	
	//output the table, k indicates the homological degree
	virtual std::set<std::pair<std::pair<int,int>,string>> output(int k, int pric){
		std::set<std::pair<std::pair<int,int>,string>> lst;
		for(auto tm : *this){
			if((tm.tag.second==Invalid || diff_length(tm)!=0) && filtration(tm.cycle)<pric)
				lst.emplace(output(tm,k));
		}
		return lst;
	}
	
	//construct a full tag from the leading term
	virtual typename SS_entry<cycle_name, ring>::value_type get_tag(cycle_name)=0;

	//construct the cycle pot. pot has additional leading number for sorting
	//reasons... When rows is non-null, only candidates built from those
	//Ptag row indices are considered (used to restrict to one degree bucket)
	virtual std::set<cycle_name> cycle_pot(int pric, const std::vector<unsigned>* rows=nullptr)=0;

	//filter the pot to delete those boundaries
	void filter_pot(std::set<cycle_name>  &pot, SS_table<cycle_name,ring> &prev_table){
		for(auto &et : prev_table){
		if(tagged(et))
			pot.erase(et.cycle);
		}
	}

	//same as filter_pot, but iterates the (small, per-bucket) pot against a
	//precomputed set of prev_table's tagged cycle names instead of iterating
	//prev_table itself -- prev_table can be large and doesn't shrink per
	//bucket, so scanning it once per bucket (as a naive per-bucket filter_pot
	//call would) is an O(#buckets x |prev_table|) regression; build the
	//tagged-name set once per homological degree instead (see make_table)
	static void filter_pot_against(std::set<cycle_name> &pot, const std::set<cycle_name> &tagged_cycles){
		for(auto it = pot.begin(); it != pot.end();){
			if(tagged_cycles.count(*it)) it = pot.erase(it);
			else ++it;
		}
	}

	//group this table's Ptag row indices by internal (t) degree, so
	//make_table can process one homological degree's search one degree
	//bucket at a time (see algNov.cpp's bucket_by_degree, which is what
	//actually implements this for algNov_table/Boc_table)
	virtual std::map<int,std::vector<unsigned>> bucket_rows()=0;

	//a fresh, empty instance of this table's own concrete type (same
	//dynamic type as *this, so filtration/naming/v_valuation etc. still
	//virtual-dispatch correctly), configured to search/build against the
	//same Ptag/Pcyc as *this. Used by make_table to give each degree
	//bucket its own small accumulation structure instead of searching the
	//whole (much larger) shared table -- see NOTES.md's bucketing entry.
	virtual SS_table<cycle_name,ring>* make_bucket_table()=0;

	//construct the table form a set of candidate of tags. Partitioned by
	//internal (t) degree: the differential is degree-preserving, so a
	//candidate can only ever cycle_index.find()-match another entry of the
	//same degree -- different buckets can never interact except through T
	//(the previous table), and different buckets' T.insert() calls always
	//touch disjoint keys, so processing buckets one at a time (in any
	//order) is equivalent to the old single flat search, just against a
	//much smaller table per candidate.
	//
	//Cross-bucket parallelism: buckets run inside one persistent
	//`#pragma omp parallel` region for this whole homological degree, not a
	//fresh thread team per bucket (that fork/join-per-call pattern was
	//already tried elsewhere in this codebase and reverted -- see NOTES.md).
	//Each thread is handed a *static* subset of buckets, decided once up
	//front by an LPT (longest-processing-time-first) partition -- sort
	//buckets descending by candidate (pot) count, greedily assign each to
	//the currently lightest thread -- since bucket sizes are skewed (up to
	//~3x at degree-270 scale) and a plain round-robin/equal-count split
	//would not balance real work.
	//
	//Two shared containers used to be mutated directly inside the
	//sequential per-bucket loop: T (via T.insert() on the Boundary branch)
	//and *this (via the temp->this merge at the end of each bucket).
	//insert() does a std::vector::push_back plus two std::map::emplace()s
	//on a container shared by every bucket -- a real data race once buckets
	//run concurrently, not a keyed operation different buckets could safely
	//interleave. Both are staged per-thread during the parallel region and
	//merged into T/*this serially afterward -- the same temp->this merge
	//pattern the sequential version already used, just pushed out one more
	//level.
	void make_table(int pric, matrix<ring>* M, SS_table& T){
		auto buckets = bucket_rows();
		//build the set of T's tagged cycle names once per homological degree
		//(not once per bucket -- see filter_pot_against) so per-bucket
		//filtering stays cheap regardless of how many buckets there are
		std::set<cycle_name> tagged_cycles;
		for(auto &et : T)
			if(tagged(et))
				tagged_cycles.insert(et.cycle);

		//build every bucket's (already-filtered) candidate pot up front --
		//needed both to drive the parallel loop below and to size the LPT
		//load-balancing partition
		std::vector<std::set<cycle_name>> pots;
		pots.reserve(buckets.size());
		for(auto &kv : buckets){
			auto pot = cycle_pot(pric, &kv.second);
			filter_pot_against(pot, tagged_cycles);
			pots.push_back(std::move(pot));
		}

		int nb = (int)pots.size();
		int nthreads = std::max(1, std::min(omp_get_max_threads(), std::max(nb,1)));

		//LPT partition: sort bucket indices descending by candidate count,
		//then greedily assign each to whichever thread currently has the
		//least assigned work -- a static decision made once, here, not a
		//dynamic "grab a bucket when idle" scheme
		std::vector<int> order(nb);
		for(int i=0;i<nb;++i) order[i]=i;
		std::sort(order.begin(), order.end(), [&](int a, int b){
			return pots[a].size() > pots[b].size(); });
		std::vector<std::vector<int>> assign(nthreads);
		std::vector<size_t> load(nthreads, 0);
		for(int idx : order){
			int t = (int)(std::min_element(load.begin(), load.end()) - load.begin());
			assign[t].push_back(idx);
			load[t] += pots[idx].size();
		}

		//per-thread staging for the two deferred shared-container merges
		//(see comment above); the only place T/*this are actually touched
		//is the serial merge loop after the parallel region closes
		std::vector<std::vector<SS_entry<cycle_name,ring>>> T_staging(nthreads), this_staging(nthreads);

		#pragma omp parallel num_threads(nthreads)
		{
			int tid = omp_get_thread_num();
			for(int bidx : assign[tid]){
				auto &pot = pots[bidx];

				//this bucket's own fresh, empty table -- entries found while
				//processing this bucket accumulate here, not in the shared *this
				std::unique_ptr<SS_table<cycle_name,ring>> temp(make_bucket_table());

				//construct the Curtis table. Note we start with tags of highest filtration
				for(auto it = pot.rbegin(); it != pot.rend(); it++){
					//get the canditate of the cycle
					auto vs = get_tag(*it);
					//a candidate whose tag vector is genuinely zero (e.g. a
					//manufactured v0^j-multiple candidate with j beyond the
					//ring's precision, so its coefficient reduced to 0) is
					//the zero module element -- not a cycle, not a boundary,
					//not any class at all. Recording it (as either branch
					//below would) manufactures a fake surviving generator
					//labeled with the raw, unbounded candidate name.
					if(Modop->isZero(vs)) continue;
					//get the boundary of the canditate
					auto bs = M->maps_to(vs);

					//try to add lower terms to make it a cycle -- searched
					//against this bucket's own (small) table, not the shared one
					auto BT = temp->simplify(bs,vs,pric);

					//make the new entry
					SS_entry<cycle_name, ring> neE;
					//if the candidate can be made a cycle, stage it for the
					//previous table T (deferred -- see comment above)
					if(BT==Boundary){
						neE = {invalid(), *it, Modop->zero(), vs};
						T_staging[tid].push_back(std::move(neE));
					} //if we find a nontrivial tag, then insert into this bucket's table
					else{
						//find the leading term if the boundary
						unsigned led = leading_term(bs);
						//make the new entry
						auto bn = naming(bs.dataArray[led].coeficient, bs.dataArray[led].ind);
						neE = {*it, bn, vs, bs};
						temp->insert(neE);
					}
				}
				//stage this bucket's entries for the real, shared table
				//(deferred -- see comment above)
				for(auto &e : *temp)
					this_staging[tid].push_back(e);
			}
		}

		//serial merge -- the only place T/*this are mutated for this call
		for(auto &v : T_staging)
			for(auto &e : v)
				T.insert(e);
		for(auto &v : this_staging)
			for(auto &e : v)
				this->insert(e);
	}
	
	//Return the set of cycles modulo filtration pric boundaries and cycles, original vector should be discarded after
	std::tuple<std::vector<cycle_name>, std::vector<cycle_name>, std::vector<cycle_name>> name_of_cycle(typename SS_entry<cycle_name, ring>::value_type& v, SS_table *next_table, int pric, bool check_cycle=true){
		//the set of tags, boundaries and cycles in v
		std::vector<cycle_name> tgs, bdrs, cycs;
		while(!Modop->isZero(v)){
			//find the leading term
			unsigned led = leading_term(v);
			cycle_name leading_name = naming(v.dataArray[led].coeficient, v.dataArray[led].ind);
			//if arrive at the desired pricision then return
			if(filtration(leading_name)>=pric) break;
			
			//find the entry for the leading term,
			auto it = cycle_index.find(leading_name);
			
			if(it == cycle_index.end()){
				//when the leading term is not a cycle, try it with the next table
				if(next_table == NULL){
					//nothing left to consult and the leading term isn't a
					//cycle here -- this used to just print and fall through,
					//leaving `v` unchanged, so the enclosing while(!isZero(v))
					//loop re-hit the exact same leading term forever: an
					//infinite loop whose only visible symptom was unbounded
					//"invalid table!:..." spam to stderr. Not a transient/
					//retryable condition -- the table this search depends on
					//is genuinely missing this candidate -- so abort loudly
					//instead of spinning.
					std::cerr << "name_of_cycle: leading term gen_ind=" << leading_name.second
						<< " (filtration=" << filtration(leading_name) << ", pric="
						<< pric << ") is not a cycle and there is no next table "
						"to check as a tag -- the table this search depends on "
						"is incomplete for this candidate.\n" << std::flush;
					std::abort();
				}
				else{
					//try to find the leading term as a tag
					auto ot = next_table->tag_index.find(leading_name);
					if(ot == next_table->tag_index.end()){
						//same non-terminating fall-through as the NULL-table
						//branch above, just one table further out -- see that
						//branch's comment.
						std::cerr << "name_of_cycle: leading term gen_ind=" << leading_name.second
							<< " (filtration=" << filtration(leading_name) << ", pric="
							<< pric << ") is not a cycle here and not a tag in "
							"the next table either -- the table this search "
							"depends on is incomplete for this candidate.\n"
							<< std::flush;
						std::abort();
					}
					else{
						if(filtration(next_table->at(ot->second).cycle)<pric && check_cycle)
							std::cerr << "invalid table!";
						//subtract r copies of the full tag, r being the
						//leading coefficient's residue mod p -- NOT the whole
						//coefficient, see cancel_residue
						auto cr = next_table->cancel_residue(v.dataArray[led].ind, v.dataArray[led].coeficient, next_table->at(ot->second).full_tag);
						v = Modop->add(v, Modop->scalor_mult(cr.first, next_table->at(ot->second).full_tag));
						//record the tag once per unit of its F_p coefficient
						for(unsigned i=0;i<cr.second;++i)
							tgs.push_back(leading_name);
					}
				}
			}
			else {
				//when found in the cycle list, subtract r copies, r being the
				//leading coefficient's residue mod p -- NOT the whole
				//coefficient, see cancel_residue
				auto cr = cancel_residue(v.dataArray[led].ind, v.dataArray[led].coeficient, this->at(it->second).full_cycle);
				v = Modop->add(v, Modop->scalor_mult(cr.first, this->at(it->second).full_cycle));
				//if it is untagged, mark as cycle, otherwise mark as boundary,
				//once per unit of its F_p coefficient
				for(unsigned i=0;i<cr.second;++i){
					if(!tagged(this->at(it->second)))
						cycs.push_back(leading_name);
					else
						bdrs.push_back(leading_name);
				}
			}
		}
		return std::make_tuple(tgs,bdrs,cycs);
	}
    
    //combine the tag and cycles
    std::pair<std::vector<cycle_name>, std::vector<cycle_name>> combine_cycles(std::tuple<std::vector<cycle_name>, std::vector<cycle_name>, std::vector<cycle_name>> tc){
		std::vector<cycle_name> cycs = std::get<0>(tc);
		for(auto tm : std::get<2>(tc))
			cycs.push_back(tm);
		return std::make_pair(std::get<1>(tc), cycs);
	}
    
	//Return the set of boundaries and cycles, original vector should be discarded after
	std::pair<std::vector<cycle_name>, std::vector<cycle_name>> name_of_cycle(typename SS_entry<cycle_name, ring>::value_type& v, int pric=64){
		auto cycls = name_of_cycle(v,NULL,pric);
		return std::make_pair(std::get<1>(cycls), std::get<2>(cycls));
	}
};
