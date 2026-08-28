//SS.h
#pragma once
#include"matrices.h"
#include<memory>

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

	//simplify a cycle by subtraction boundaries using the data in the table
	cycle_type simplify(typename SS_entry<cycle_name, ring>::value_type& cyc, typename SS_entry<cycle_name, ring>::value_type& htpy, int pric, int depth=0){
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
		return simplify(cyc,htpy,pric,depth+1);
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
	void make_table(int pric, matrix<ring>* M, SS_table& T){
		auto buckets = bucket_rows();
		//build the set of T's tagged cycle names once per homological degree
		//(not once per bucket -- see filter_pot_against) so per-bucket
		//filtering stays cheap regardless of how many buckets there are
		std::set<cycle_name> tagged_cycles;
		for(auto &et : T)
			if(tagged(et))
				tagged_cycles.insert(et.cycle);
		for(auto &kv : buckets){
			auto &rows = kv.second;
			//construct the cycle_pot restricted to this bucket's rows
			auto pot = cycle_pot(pric, &rows);
			//filter the cycle_pot against T's (precomputed) tagged names
			filter_pot_against(pot, tagged_cycles);

			//this bucket's own fresh, empty table -- entries found while
			//processing this bucket accumulate here, not in the shared *this
			std::unique_ptr<SS_table<cycle_name,ring>> temp(make_bucket_table());

			//construct the Curtis table. Note we start with tags of highest filtration
			for(auto it = pot.rbegin(); it != pot.rend(); it++){
				//get the canditate of the cycle
				auto vs = get_tag(*it);
				//get the boundary of the canditate
				auto bs = M->maps_to(vs);

				//try to add lower terms to make it a cycle -- searched
				//against this bucket's own (small) table, not the shared one
				auto BT = temp->simplify(bs,vs,pric);

				//make the new entry
				SS_entry<cycle_name, ring> neE;
				//if the candidate can be made a cycle, then insert into the previous table
				if(BT==Boundary){
					neE = {invalid(), *it, Modop->zero(), vs};
					T.insert(neE);
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
			//merge this bucket's entries into the real, shared table
			for(auto &e : *temp)
				this->insert(e);
		}
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
				if(next_table == NULL)
					std::cerr << "invalid table!:end of table";
				else{
					//try to find the leading term as a tag
					auto ot = next_table->tag_index.find(leading_name);
					if(ot == next_table->tag_index.end())
						std::cerr << "invalid table!: no such tag";
					else{
						if(filtration(next_table->at(ot->second).cycle)<pric && check_cycle)
							std::cerr << "invalid table!";
						//subtract the full tag, scaled to cancel the leading term exactly (see cancel_scalar)
						auto scalar = next_table->cancel_scalar(v.dataArray[led].ind, v.dataArray[led].coeficient, next_table->at(ot->second).full_tag);
						v = Modop->add(v, Modop->scalor_mult(scalar, next_table->at(ot->second).full_tag));
						//add the information of the tag
						tgs.push_back(leading_name);
					}
				}
			}
			else {
				//when find it in the cycle list, modify the vector, scaled to cancel the leading term exactly (see cancel_scalar)
				auto scalar = cancel_scalar(v.dataArray[led].ind, v.dataArray[led].coeficient, this->at(it->second).full_cycle);
				v = Modop->add(v, Modop->scalor_mult(scalar, this->at(it->second).full_cycle));
				//if it is untagged, mark as cycle, otherwise mark as boundary
				if(!tagged(this->at(it->second))) 
					cycs.push_back(leading_name);
				else 
					bdrs.push_back(leading_name);
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
