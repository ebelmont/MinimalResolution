//multiplication.cpp
#include"multiplication.h"

//make the multiplication table for the primitives
void multiplication::make_eta_R_multiplier(BPBP const& x,matrix<BP>* result, int deg){
	std::cout << "making multiplier!" << std::flush;
	std::function<vectors<matrix_index, BP>(int)> rows = [this,&x,deg](int i){
		//filter out those elements beyond the scope of the data
		if(BPoper->mon_index.max_degree < deg or i > (int)BPoper->mon_index.ranksBelow[BPoper->mon_index.max_degree-deg])
			return BPoper->BPMod_opers.zero();
		//construct v^e
		auto v0 = BPoper->singleton(BPoper->mon_index.mon_array[i]);
		//the primitive is etaR(v^e)
		auto a0 = BPoper->etaR(v0);
		//compute etaR(v^e)x
		auto m0 = BPoper->BPBP_opers.multiply(a0,x);
		//transform to a vector
		return BPoper->algebroid2vector(m0,0);
	};
	//construct the matrix
	result->construct(BPoper->mon_index.number_of_all_mons(),rows);
}

//construct the multiplication matrix
void multiplication::multily_matrix(int max_deg, matrix<BP> *multiplier, matrix<BP> *res_map, FreeBPCoMod* F, primitive_data* prim_next, matrix<Z3> *result){
	std::function<vectors<matrix_index,Z3>(int)> rows = [max_deg, this, multiplier, res_map, F, prim_next](int i){
		//filter out those out of range
		if(F->degree(i)>max_deg) return BPoper->Z3Mod_opers.zero();
		//compute d(etaR(v^e)x)
		auto v = res_map->maps_to(F->multiply_using_table(i,*multiplier));
		return prim_next->expand(v);
	};
	result->construct(F->rank(),rows);
}

//make multiplication table
template<typename cycle_name, typename ring>
multiplication_table<cycle_name> multiplication::mult_extension(matrix<ring> *mult_matrix, SS_table<cycle_name,ring> &cur_table, SS_table<cycle_name,ring> &next_table, SS_table<cycle_name,ring> &third_table, int pric){
	multiplication_table<cycle_name> res;
	//compute the extensions for the entries in the table
	for(auto &tm : cur_table){
		//skip invalid entries
		if(!cur_table.valid(tm.cycle)) continue;
		//skip the tagged entries
		if(cur_table.tagged(tm)) continue;
		//get the image of the full cycle under multiplication
		auto v = mult_matrix->maps_to(tm.full_cycle);
		//construct the new entry
		auto cycls = next_table.name_of_cycle(v,&third_table,pric);
		auto cm = next_table.combine_cycles(cycls);
		multiplication_table_entry<cycle_name> nm = {tm.cycle,cm};
		res.push_back(nm);
	}
	
	//then deal with the tags
	for(auto &tm : next_table){
		//skip invalid entries
		if(!next_table.valid(tm.cycle)) continue;
		//skip untagged entries
		if(!next_table.tagged(tm)) continue;
		//skip tags with non-trivial boundaries
		if(next_table.filtration(tm.cycle)<=pric) continue;
		auto v = mult_matrix->maps_to(tm.full_tag);
		//construct new entry
		auto cycls = next_table.name_of_cycle(v,&third_table,pric);
		auto cm = next_table.combine_cycles(cycls);
		multiplication_table_entry<cycle_name> nm = {tm.tag,cm};
		res.push_back(nm);
	}
	return res;
}

//make multiplication table
template<typename cycle_name, typename ring>
multiplication_table<cycle_name> multiplication::mult_extension1(matrix<ring> *mult_matrix, SS_table<cycle_name,ring> &cur_table, SS_table<cycle_name,ring> &next_table, SS_table<cycle_name,ring> &third_table, int pric){
	multiplication_table<cycle_name> res;
	int counter = 0;
	//compute the extensions for the entries in the table
	for(auto &tm : cur_table){
		std::cout << "\r" << counter++;
		//skip invalid entries
		if(!cur_table.valid(tm.cycle)) continue;
		//skip the tagged entries
		if(cur_table.tagged(tm)) continue;
		//skip
		if(cur_table.filtration(tm.cycle)>pric) continue;
		//get the image of the full cycle under multiplication
		auto v = mult_matrix->maps_to(tm.full_cycle);
		//construct the new entry
		auto cycls = next_table.name_of_cycle(v,&third_table,pric);
		auto cm = next_table.combine_cycles(cycls);
		multiplication_table_entry<cycle_name> nm = {tm.cycle,cm};
		res.push_back(nm);
	}
	
	//then deal with the tags
	for(auto &tm : next_table){
		std::cout << "\r" << counter++ << "/" << next_table.size();
		//skip invalid entries
		if(!next_table.valid(tm.cycle)) continue;
		//skip untagged entries
		if(!next_table.tagged(tm)) continue;
		//skip
		std::cout << cur_table.filtration(tm.tag);
		if(cur_table.filtration(tm.tag)>pric) continue;
		auto v = mult_matrix->maps_to(tm.full_tag);
		//construct new entry
		auto cycls = next_table.name_of_cycle(v,&third_table,pric,false);
		auto cm = next_table.combine_cycles(cycls);
		multiplication_table_entry<cycle_name> nm = {tm.tag,cm};
		res.push_back(nm);
	}
	return res;
}

//make multiplacation table from the algebraic Novikov table of a complex
std::vector<multiplication_table<cycle_name>> multiplication::mult_extension(matrix<BP> *multiplier, int max_deg, int resolution_length, std::vector<FreeBPCoMod>& generators, string maps_filename, BPComplex& Cm, algNov_tables Tb, matrix<BP> *mp, matrix<Z3> *mm, int pric, bool fixedpric){
	//open the files of the maps
	std::fstream maps_file(maps_filename, std::ios::in | std::ios::binary);
	std::vector<multiplication_table<cycle_name>> result;
	//make the tables for each term in the resolution	
	for(int i=1;i<resolution_length-1;++i){
		//load the maps in the resolution
		mp->load(maps_file);
		//compute the matrix for the multiplacation
		mm->clear();
		multily_matrix(max_deg, multiplier, mp, &generators[i-1], &Cm.Prims[i], mm);
		//compute the multiplication table for the current term in the resolution.
		//no filtration cap at any homological degree (see BP_init.cpp's
		//NO_FILTRATION_LIMIT); it used to be L - s, capping s+f.
		int cpric = pric;
		multiplication_table<cycle_name> new_table = mult_extension<cycle_name,Z3>(mm, *(Tb.tables)[i-1], *(Tb.tables)[i], *(Tb.tables)[i+1], cpric);
		result.push_back(new_table);
	}
	return result;
}

//make multiplacation table from the algebraic Novikov table of a complex
std::vector<multiplication_table<cycle_name>> multiplication::mult_extension1(matrix<BP> *multiplier, int max_deg, int resolution_length, std::vector<FreeBPCoMod>& generators, string maps_filename, BPComplex& Cm, algNov_tables Tb, matrix<BP> *mp, matrix<Z3> *mm, int pric, bool fixedpric){
	//open the files of the maps
	std::fstream maps_file(maps_filename, std::ios::in | std::ios::binary);
	std::vector<multiplication_table<cycle_name>> result;
	//make the tables for each term in the resolution	
	for(int i=1;i<resolution_length-1;++i){
		std::cout << "\n" << i << std::flush;
		//load the maps in the resolution
		mp->load(maps_file);
		//compute the matrix for the multiplacation
		mm->clear();
		multily_matrix(max_deg, multiplier, mp, &generators[i-1], &Cm.Prims[i], mm);
		//compute the multiplication table for the current term in the resolution.
		//no filtration cap at any homological degree (see BP_init.cpp's
		//NO_FILTRATION_LIMIT); it used to be L - s, capping s+f.
		int cpric = pric;
		multiplication_table<cycle_name> new_table = mult_extension1<cycle_name,Z3>(mm, *(Tb.tables)[i-1], *(Tb.tables)[i], *(Tb.tables)[i+1], cpric);
		result.push_back(new_table);
	}
	return result;
}

//output the table
template<typename cycle_name, typename ring>
string multiplication::output_multiplication_table(multiplication_table<cycle_name> const &MT, int k, int shift, int pric, SS_table<cycle_name, ring> &Tb){
	string res;
	for(auto &tm: MT){
		//skip those entries out of range
		if(Tb.filtration(tm.original_name) >= pric) continue;
		//output the origin
		res+= Tb.output(tm.original_name,k) + "\t->\t";
		//output the result
		for(auto sm: tm.multiplied_names.second)
			res+= Tb.output(sm,k+shift)+"+";
		res+= "o\n";
	}
	return res;
}

//output the tables
string multiplication::output_multiplication_table(std::vector<multiplication_table<cycle_name>> const &MTs, int shift, int pric){
	algNov_table Tb;
	string res;
	for(unsigned i=0;i<MTs.size();++i)
		res+= output_multiplication_table(MTs[i],i,shift, pric,Tb);
	return res;
}

//output the table
template<typename cycle_name, typename ring>
string multiplication::output_multiplication_table1(multiplication_table<cycle_name> const &MT, int k, int shift, int pric, SS_table<cycle_name, ring> &Tb){
	string res;
	for(auto &tm: MT){
		//skip those entries out of range
		if(Tb.filtration(tm.original_name) >= pric) continue;
		//output the origin
		res+= Tb.output(tm.original_name,k) + "\t->\t";
		//output the result
		for(auto sm: tm.multiplied_names.second)
			res+= Tb.output(sm,k+shift)+"+";
		for(auto sm: tm.multiplied_names.first)
			res+= Tb.output(sm,k+shift)+"+";
		res+= "o\n";
	}
	return res;
}

//output the tables
string multiplication::output_multiplication_table1(std::vector<multiplication_table<cycle_name>> const &MTs, int shift, int pric){
	algNov_table Tb;
	string res;
	for(unsigned i=0;i<MTs.size();++i)
		res+= output_multiplication_table1(MTs[i],i,shift, pric,Tb);
	return res;
}

//the extension by 3
template<typename cycle_name, typename ring>
multiplication_table<cycle_name> multiplication::three_extension(SS_table<cycle_name,ring>& cur_table, int pric){
	multiplication_table<cycle_name> res;
	//the element 3
	static auto three = cur_table.Modop->ringOper->unit(3);
	//compute the 3-extensions
	for(auto &tm : cur_table){
		//skip invalid entries
		if(!cur_table.valid(tm.cycle)) continue;
		//Used to also skip tagged (differential-target) entries here. For
		//algNov_table that's harmless (only permanent-cycle names are ever
		//looked up downstream, since deg_table/BPAANSS_table.txt only lists
		//permanent cycles), but for Boc_table it silently left every
		//differential target with no "multiply by 3" entry in
		//BPBocSS_a0.txt at all -- delta() needs exactly that (the target of
		//a Bockstein differential, y in `y <- x | dr`, is itself frequently
		//v0-divisible and needs dividing before B2A can translate it).
		//Confirmed by build+run test (degree 20 repro): removing this skip
		//adds only new a0 entries -- every entry present before is
		//byte-identical after, both for algNov_table's own a0.txt and for
		//Boc_table's. See NOTES.md.
		//compute the multiplication of the cycle by 3
		auto v = cur_table.Modop->scalor_mult(three, tm.full_cycle);
		//compute the name of the result
		multiplication_table_entry<cycle_name> nm = {tm.cycle,cur_table.name_of_cycle(v,pric)};
		res.push_back(nm);
	}
	return res;
}

//make 3-extension table from the algebraic Novikov table of a complex
std::vector<multiplication_table<cycle_name>> multiplication::three_extension(int resolution_length, BPComplex& Cm, algNov_tables Tb, int pric){
	std::vector<multiplication_table<cycle_name>> result;
	//make the tables for each term in the resolution	
	for(int i=0;i<resolution_length;++i){
		//compute the multiplication table for the current term in the resolution.
		//no filtration cap at any homological degree (see BP_init.cpp's
		//NO_FILTRATION_LIMIT); it used to be L - s, capping s+f.
		multiplication_table<cycle_name> new_table = three_extension(*(Tb.tables)[i], pric);
		result.push_back(new_table);
	}
	return result;
}

//constructor
multiplication::multiplication(BP_Op* BP_oper){
	BPoper = BP_oper; }
