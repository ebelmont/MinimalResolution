//Boc.cpp
#include"Boc.h"
#include"fast_itos.h"

int Boc_table::v_valuation(Z3 r, exponent e){
	return v_valuation_from(Z3Oper->valuation(r), e);
}

//see algNov.h/algNov.cpp -- same, given a precomputed valuation
int Boc_table::v_valuation_from(unsigned val, exponent){
	return (int)val;
}

void Boc_tables::set_table(std::vector<algNov_table>*){
	std::cerr << "action not allowed!\n";
}

void Boc_tables::set_table(std::vector<Boc_table>* tbs){
	tables.resize(tbs->size());
	Btables.resize(tbs->size());
	for(unsigned i=0; i<tbs->size(); ++i){
		tables[i] = &tbs->at(i);
		Btables[i] = &tbs->at(i);
	}
}

int Boc_table::num_v(cycle_name cyc){
	int res = 0;
	for(unsigned i=1; i<cyc.first.size(); ++i)
		res += cyc.first[i];
	return res;
}

std::pair<std::pair<int,int>,string> Boc_table::output(SS_entry<cycle_name, Z3>& et, int k){
	string res;
	//if it is a untagged cycle
	if(!tagged(et)) 
		res = this->algNov_table::output(et.cycle,k);
	//if it is tagged, output cyc <- tag
	else 
		res = this->algNov_table::output(et.cycle,k) + "\t<-\t" + this->algNov_table::output(et.tag,k-1) + "\t|d" + fast_itos(diff_length(et));
	
	std::pair<int,int> degs = std::make_pair(degree(et.cycle)*2-k, k+filtration(et.cycle));
	return std::make_pair(degs, res + "\t|deg=" + fast_itos(degs.first));
}

std::set<std::pair<std::pair<int,int>,string>> Boc_table::output_table(int k, int pric){
	std::set<std::pair<std::pair<int,int>,string>> lst;
	for(auto tm : *this){
		//Used to also skip any entry with nonzero v0-content here ("skip
		//those entries with v0"). That silently dropped real Bockstein-space
		//classes that 185_BPBocSS_theta*.txt/_h0.txt reference directly (as
		//sources and as target terms) -- delta() (beta.py/convert_e5.py)
		//needs exactly these as lookups in 185_BPBocSS_table.txt, and with
		//the filter in place they were simply absent, not derivable from any
		//other file. Confirmed by direct build+run test (degree 20 repro)
		//that removing the filter recovers them with zero change to any
		//previously-printed (v0-free) entry. See NOTES.md. Now that the theta
		//tables are named through v0^0 only (mult_theta), delta() no longer
		//needs these rows, but convert_e5.py's first-copy/second-copy
		//classification does: second-copy classes are v0-divisible survivors
		//(v0^{40-k} x), and with the filter back it finds none of them.
		if(filtration(tm.cycle)<pric){
			lst.emplace(output(tm,k));
		}
	}
	return lst;
}


//transform Bockstein name into Algebraic Novikov name
multiplication_table<cycle_name> Boc_table::Bname2Aname(algNov_table &Atable, int pric, bool simple){
	//initialize the resulting table
	multiplication_table<cycle_name> res;
	//go through all the entries
	for(auto &tm : *this){
		//skip those tagged entries
		if(tagged(tm)) continue;
		//skip those entries divizible by v0 if we just want simple tables
		if(simple && tm.cycle.first[0]>0) continue;
		//read the full cycle
		auto v = tm.full_cycle;
		//construct the new entry
		multiplication_table_entry<cycle_name> nm = {tm.cycle, Atable.name_of_cycle(v, pric)};
		res.push_back(nm);
	}
	return res;
}

//transform Bockstein names into Algebraic Novikov names
std::vector<multiplication_table<cycle_name>> Boc_tables::Bname2Anames(int resolution_length, algNov_tables& AT, int pric){
	std::vector<multiplication_table<cycle_name>> result;
	//see algNov.cpp::table_of_complex's identical fix for why min(...)
	//replaces the old pric-i (pric is now the ring-precision ceiling, not
	//resolution_length, so decrementing it by i re-shrinks an
	//already-clamped bound by an unrelated quantity).
	for(int i=0; i<resolution_length; ++i){
		result.push_back(Btables[i]->Bname2Aname(*(AT.tables[i]), std::min(resolution_length-i, pric)));
	}
	return result;   
}

//constructor
Boc_table::Boc_table() : SS_table(&algNov_table::BPoper->Z3Mod_opers){}

//see SS.h::make_table / algNov.cpp::algNov_table::make_bucket_table -- must
//return a genuine Boc_table (not algNov_table) so v_valuation/
//v_valuation_from keep virtual-dispatching to the Bockstein override
SS_table<cycle_name,Z3>* Boc_table::make_bucket_table(){
	Boc_table *t = new Boc_table();
	t->Ptag = Ptag;
	t->Pcyc = Pcyc;
	return t;
}
