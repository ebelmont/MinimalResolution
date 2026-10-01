//algNov.cpp
#include"algNov.h"
#include"fast_itos.h"
#include<cstdio>

//the operators
Z3_Op *algNov_table::Z3Oper;
ModuleOp<matrix_index, Z3> *algNov_table::Modop;
BP_Op *algNov_table::BPoper;

//filtration of a cycle
int algNov_table::filtration(cycle_name cyc){
	return cyc.first[0];
}

//the filtration by number of v's
int algNov_table::v_valuation(Z3 r, exponent e){
	return v_valuation_from(Z3Oper->valuation(r), e);
}

//same as v_valuation, but given an already-computed Z3Oper->valuation(r),
//so callers that separately need the raw valuation (naming(),
//leading_term()) can avoid computing it twice -- Z3_Op::valuation averages
//~19 divide-steps/call on degree-270 data, so this matters on the
//leading_term hot path. Inlined rather than routed through the generic
//total_deg(e, [](int){return 1;}) to avoid a std::function indirection per call.
int algNov_table::v_valuation_from(unsigned val, exponent e){
	int s = 0;
	for(int i=1; i<=maxVar; ++i)
		s += xnVal(e,i);
	return (int)val + s;
}
	
//transform a term into a cycle name
cycle_name algNov_table::naming(const Z3 &coeficient, matrix_index ind, primitive_data &P){
	//get the exponent
	exponent e = P[ind].coeficient;
	//construct the cycle
	std::vector<int> eo;
	//final size is always exactly maxVar+2, known up front
	eo.reserve(maxVar+2);
	//compute the valuation once and reuse it for both the algebraic Novikov
	//filtration and the raw v0-count, instead of calling
	//Z3Oper->valuation(coeficient) twice (see v_valuation_from)
	unsigned val = Z3Oper->valuation(coeficient);
	//first compute the algebraic Novikov filtration
	eo.push_back(v_valuation_from(val,e));
	//then get the number of v0's
	eo.push_back(val);
	//then the exponent
	for(int i=1; i<=maxVar; ++i)
		eo.push_back(xnVal(e,i));
	
	//the cycle
	return {eo, P[ind].gen_pos};
}

//transform a term into a cycle name, using Pcyc
cycle_name algNov_table::naming(const Z3 &coeficient, matrix_index ind){
	return naming(coeficient, ind, *Pcyc);
}
	
//computing the leading term of a cycle vector
//compares terms by a fixed-size stack key (no heap allocation) instead of
//calling naming() -- a heap-allocating pair<vector<int>,int> -- per
//candidate, and caches the current-best key rather than recomputing it
//each iteration. The key's fields and order exactly match naming()'s eo
//layout (v_valuation, v0-valuation, xnVal(1..maxVar)) plus gen_pos as the
//tiebreaker (naming()'s pair comparison falls back to .second on a tie in
//.first), so the winning index is unchanged. naming() itself is still
//called exactly once per leading_term() call, by simplify() right after
//this returns, to materialize the winner's actual name.
unsigned algNov_table::leading_term(SS_entry<cycle_name, Z3>::value_type& v){
	//case of zero
	if(v.size()==0) return -1;

	struct Key{ int f[maxVar+2]; int gen_pos; };
	auto compute_key = [&v,this](unsigned idx){
		matrix_index ind = v.dataArray[idx].ind;
		auto &pe = (*Pcyc)[ind];
		exponent e = pe.coeficient;
		unsigned val = Z3Oper->valuation(v.dataArray[idx].coeficient);
		Key k;
		k.f[0] = v_valuation_from(val,e);
		k.f[1] = (int)val;
		for(int i=1; i<=maxVar; ++i)
			k.f[i+1] = xnVal(e,i);
		k.gen_pos = pe.gen_pos;
		return k;
	};
	auto key_less = [](const Key &a, const Key &b){
		for(int i=0; i<maxVar+2; ++i)
			if(a.f[i]!=b.f[i]) return a.f[i]<b.f[i];
		return a.gen_pos < b.gen_pos;
	};

	//find the leading term with the smallest cycle name
	unsigned led = 0;
	Key best = compute_key(0);
	for(unsigned i=1; i<v.size(); ++i){
		Key ki = compute_key(i);
		if(key_less(ki,best)){
			led = i;
			best = ki;
		}
	}
	return led;
}
	
//set the operations
void algNov_table::set_op(BP_Op *BP_oper){
	BPoper = BP_oper;
	Modop = &BPoper->Z3Mod_opers;
	Z3Oper = BPoper->Z3_oper;
}
	
//constructors
algNov_table::algNov_table() : SS_table(&BPoper->Z3Mod_opers){}

//check the entry is a tag or not
bool algNov_table::tagged(const SS_entry<cycle_name,Z3> &et){
	return et.tag.second != Invalid;
}

cycle_name algNov_table::invalid(){
	std::vector<int> vd;
	return std::make_pair(vd,Invalid);
}

//construct the cycle pot. pot has additional leading number for sorting
//reasons... When rows is given, only those Ptag row indices are considered
//(used to restrict the search to one internal-degree bucket -- see
//bucket_rows()/SS.h::make_table); passing nullptr reproduces the original
//whole-Ptag behavior.
std::set<cycle_name> algNov_table::cycle_pot(int pric, const std::vector<unsigned>* rows){
	std::set<cycle_name> res;
	auto add_row = [&](unsigned i){
		//construct ((l,0,v1^e1,...),gen_ind)
		auto s = naming(Z3Oper->unit(1), i, *Ptag);
		//compute the filtration
		int l = filtration(s);
		for(int j=0; j+l<=pric; ++j){
			//construct ((l+j, j, v1^e1, ...), gen_ind)
			s.first[0] = j+l;
			s.first[1] = j;
			res.emplace(s);
		}
	};
	if(rows){
		for(unsigned i : *rows) add_row(i);
	} else {
		for(unsigned i=0; i<Ptag->size(); ++i) add_row(i);
	}
	return res;
}

//group a primitive_data's row indices by internal (t) degree -- the
//differential is degree-preserving (v0-shifts never cross t-degree), so
//candidates built from rows of different degrees can never
//cycle_index.find()-match each other. Grouping by the existing gen_deg
//field (no new bookkeeping) is what SS.h::make_table uses to process one
//homological degree's search a bucket at a time instead of as one flat
//search over the whole table.
std::map<int,std::vector<unsigned>> bucket_by_degree(primitive_data &P){
	std::map<int,std::vector<unsigned>> buckets;
	for(unsigned i=0; i<P.size(); ++i)
		buckets[P.gen_deg[i]].push_back(i);
	return buckets;
}

std::map<int,std::vector<unsigned>> algNov_table::bucket_rows(){
	return bucket_by_degree(*Ptag);
}

//a fresh, empty algNov_table configured with the same Ptag/Pcyc as *this --
//see SS.h::make_table for why each degree bucket needs its own small,
//empty accumulation structure rather than searching the whole shared table
SS_table<cycle_name,Z3>* algNov_table::make_bucket_table(){
	algNov_table *t = new algNov_table();
	t->Ptag = Ptag;
	t->Pcyc = Pcyc;
	return t;
}

//output of cycle names
string algNov_table::output(cycle_name sn, int k){
	string res;
	for(unsigned i=1; i<sn.first.size(); ++i)
		if(sn.first[i]!=0)
			res += "v" + fast_itos(i-1) + "^" + fast_itos(sn.first[i]);
	res += "[" + fast_itos(k) + "-" + fast_itos(sn.second) + "]";
	return res;
}

//IO of entries
std::pair<std::pair<int,int>,string> algNov_table::output(SS_entry<cycle_name, Z3>& et, int k){
	string res;
	//if it is a untagged cycle
	if(!tagged(et)) 
		res = output(et.cycle,k);
	//if it is tagged, output cyc <- tag
	else 
		res = output(et.cycle,k) + "\t<-\t" + output(et.tag,k-1) +"\t|d" + fast_itos(diff_length(et)+1);
	//degree = (t-s, s+i)
	std::pair<int,int> degs = std::make_pair(degree(et.cycle)-k, k+filtration(et.cycle));
	
	return std::make_pair(degs, res + "\t|deg=(" + fast_itos(degs.first) + "," + fast_itos(degs.second)+ ")");
}

//save the cycle  name
void algNov_table::save(cycle_name nm, std::iostream& fl){
	unsigned ls = nm.first.size();
	fl.write((char*)&ls, sizeof(unsigned));
	fl.write((char*)nm.first.data(), ls*sizeof(int));
	fl.write((char*)&nm.second, sizeof(int));
}

//save the entry
void algNov_table::save(SS_entry<cycle_name, Z3>& et, std::iostream& fl){
	save(et.tag,fl);
	save(et.cycle,fl);
	Modop->save(et.full_tag,fl);
	Modop->save(et.full_cycle,fl);
}

//load the cycle name
cycle_name algNov_table::load(std::iostream& fl){
	unsigned ls;
	fl.read((char*)&ls, sizeof(unsigned));
	cycle_name res;
	res.first.resize(ls);
	fl.read((char*)res.first.data(), ls*sizeof(int));
	fl.read((char*)&res.second, sizeof(int));
	return res;
}

//load the entry
void algNov_table::load(SS_entry<cycle_name, Z3>& et, std::iostream& fl){
	et.tag = load(fl);
	et.cycle = load(fl);
	et.full_tag = Modop->load(fl);
	et.full_cycle = Modop->load(fl);
}

//construct a vector using the primitive data
SS_entry<cycle_name, Z3>::value_type algNov_table::make_vec(cycle_name cyc, primitive_data& P, ModuleOp<matrix_index,Z3> *Modop, Z3_Op *Z3Oper){ 
	//construct the exponent
	auto e = pack(cyc.first.data()+2);
	//find the entry
	auto n = P.prim_index[{(matrix_index)cyc.second, e}];
	//construct the term v0^e0v1^e1...[gen_pos]
	return Modop->singleton(n,Z3Oper->power_p(cyc.first[1]));
}

//construct a full tag from the leading term
SS_entry<cycle_name, Z3>::value_type algNov_table::get_tag(cycle_name tag){
	return make_vec(tag, *Ptag, Modop, Z3Oper); }

//the degree of cycle
int algNov_table::degree(cycle_name cyc){
	//get the exponent of cyc = ((fil, e0,e1,...),gen_pos)
	exponent e = pack(cyc.first.data()+2);
	//find the position of v^e[gen_pos]
	int n = Pcyc->prim_index[{(matrix_index)cyc.second,e}];
	//get the degre
	return Pcyc->gen_deg[n];
}
    
//construct the table using the complex
void  algNov_tables::table_of_complex(BPComplex& Comp, int pric, string back_up_prefix){
	//initialize primitive data
	tables[0]->Pcyc = &Comp.Prims[0];
	for(unsigned i=1; i<Comp.size(); ++i){
		tables[i]->Ptag = &Comp.Prims[i-1];
		tables[i]->Pcyc = &Comp.Prims[i];
	}

	//resume from whatever checkpoints already exist. table i's checkpoint
	//(written below, right after it's built) holds both tables[i-1] --
	//which is finalized at that point, since the only thing that can ever
	//add to a table after it's built is the *next* table's construction
	//(make_table's Boundary branch inserts into the previous table -- see
	//NOTES.md) -- and tables[i] itself, in its current (not yet finalized)
	//state, which is exactly what's needed to resume building from there.
	//
	//Each slot is preceded by a 1-byte flag: 1 means the slot's full
	//SS_table payload follows (as before); 0 means "unchanged since the
	//last time this table index was written" -- SS_table::insert only
	//ever appends (see NOTES.md), so an unchanged entry-count means
	//unchanged content, and the in-memory table object (already correct
	//from the previous iteration's write, or from loading it as the
	//*other* slot one file ago) can simply be left alone. `last_saved`
	//tracks, per table index, the entry-count as of its last write, so
	//both the resume-load below and the write loop further down can make
	//the same skip decision.
	std::vector<size_t> last_saved(Comp.size(), (size_t)-1);
	unsigned start = 1;
	if(back_up_prefix!=""){
		unsigned k=1;
		for(; k<Comp.size(); ++k){
			std::fstream fk(back_up_prefix + std::to_string(k), std::ios::in | std::ios::binary);
			if(!fk.is_open()) break;
			char flag;
			fk.read(&flag,1);
			if(flag){
				tables[k-1]->SS_table<cycle_name,Z3>::load(fk);
				last_saved[k-1] = tables[k-1]->size();
			}
			fk.read(&flag,1);
			if(flag){
				tables[k]->SS_table<cycle_name,Z3>::load(fk);
				last_saved[k] = tables[k]->size();
			}
		}
		if(k>1){
			std::cout << "resuming algNov tables from table " << (k-1) << "\n" << std::flush;
			start = k;
		}
	}

	//construct tables
	for(unsigned i=start; i<Comp.size(); ++i){
		std::cout << "\n" << "table:" << i << "\n" << std::flush;
		//this table's own search budget: `pric-i` reused the pre-fix
		//meaning of `pric` (== resolution_length, so pric-i was really
		//"resolution steps remaining"), but callers now pass the
		//Z3_Op::MAX_VALUATION-clamped ring-precision ceiling as `pric`
		//instead (see NOTES.md "2026-09-10: correction") -- decrementing
		//*that* by i re-shrinks an already-clamped ring bound by
		//homological degree, which is a different, unrelated quantity, and
		//wrongly excises real surviving generators once i approaches pric
		//(see NOTES.md's `three_extension`/`pric-i` entry). Comp.size()
		//is this table_of_complex's resolution_length (Complex.load_matrix
		//is called with resolution_length, giving Comp this many entries),
		//so min(Comp.size()-i, pric) reproduces the old pric-i exactly
		//when pric was still unclamped (pric==Comp.size()), and otherwise
		//correctly stays at the ring ceiling instead of shrinking with i.
		int table_pric = std::min((int)Comp.size()-(int)i, pric);

		//pre-warm the power_p cache serially, before make_table's parallel
		//region: Z3Oper is one object shared by every algNov_table/
		//Boc_table, and its power_p_cache (Z3.h) is a plain std::vector,
		//grown in place by push_back -- safe to call concurrently only
		//once every index it will ever be asked for during this table's
		//build is already present, since a v0-exponent or valuation used
		//by naming()/make_vec()/cancel_scalar() can never exceed this
		//table's own pric (the algorithm's existing filtration-vs-pric
		//invariant, not a new one introduced by parallelism). Once warmed
		//to pric, power_p(n<=pric) is a pure read of an already-sized
		//vector, safe for concurrent callers.
		algNov_table::Z3Oper->power_p(table_pric);
		tables[i]->make_table(table_pric,Comp.Maps[i],*tables[i-1]);

		//checkpoint tables[i-1] (now finalized) and tables[i] (current)
		//together in one file, written to a temp name and renamed into
		//place atomically -- so a crash mid-write leaves the last good
		//checkpoint intact instead of a half-written, silently-wrong one.
		if(back_up_prefix!=""){
			string fname = back_up_prefix + std::to_string(i);
			string tmp = fname + ".tmp";
			{
				std::fstream fk(tmp, std::ios::out | std::ios::binary);
				//skip re-serializing a slot whose entry-count hasn't
				//changed since it was last written (see comment above the
				//resume-load loop for why this is safe)
				size_t sz0 = tables[i-1]->size();
				char flag0 = (sz0 != last_saved[i-1]) ? 1 : 0;
				fk.write(&flag0,1);
				if(flag0){
					tables[i-1]->SS_table<cycle_name,Z3>::save(fk);
					last_saved[i-1] = sz0;
				}
				size_t sz1 = tables[i]->size();
				char flag1 = (sz1 != last_saved[i]) ? 1 : 0;
				fk.write(&flag1,1);
				if(flag1){
					tables[i]->SS_table<cycle_name,Z3>::save(fk);
					last_saved[i] = sz1;
				}
			}
			std::rename(tmp.c_str(), fname.c_str());
		}
	}
}

//out put the algebraic Novikov table
string algNov_tables::output_tables(){
	string res;
	std::set<std::pair<std::pair<int,int>,string>> lst;
	for(unsigned i=0; i<tables.size(); ++i){
		auto jst = tables[i]->output_table(i,tables.size()-i-1);
		for(auto sm: jst)
			lst.insert(sm);
	}
	for(auto sm : lst)
		res+=sm.second+"\n";
	return res;
}

//save the tables
void algNov_tables::save(std::iostream& fl){
	for(auto tm:tables)
		tm->SS_table<cycle_name,Z3>::save(fl);
}

//load the tables
void algNov_tables::load(std::iostream& fl){
	for(auto tm:tables)
		tm->SS_table<cycle_name,Z3>::load(fl);
}

//save to file
void algNov_tables::save(string filename){
	std::fstream fk(filename, std::ios::out | std::ios::binary);
	save(fk);
}

//load from filte
void algNov_tables::load(string filename){
	std::fstream fk(filename, std::ios::in | std::ios::binary);
	load(fk);
}
    
//set the data of tables
void algNov_tables::set_table(std::vector<algNov_table>* tbs){
	tables.resize(tbs->size());
	for(unsigned i=0; i<tbs->size(); ++i)
		tables[i] = &tbs->at(i);
}

//check if a name is valid
bool algNov_table::valid(const cycle_name& cyc){
	return cyc.second != Invalid;
}

std::set<std::pair<std::pair<int,int>,string>> algNov_table::output_table(int k, int pric){
	return SS_table<cycle_name,Z3>::output(k,pric);
}
