//trunc_init.cpp
#include"trunc_init.h"

//the initialization
TruncInit::TruncInit(int prime, int max_deg, int res_length) : trunc_oper(max_deg, prime){
	resolution_length = res_length;

	std::cout << "maximal degree:" << trunc_oper.maxDeg << "\n";
	std::cout << "resolution length:" << res_length << "\n";

	//Gamma's comultiplication is computed directly by trunc_oper.delta (a closed-form
	//binomial formula), so unlike the Steenrod algebra there is no co-product table to
	//build or load from disk.
	trunc_oper.initialize();

	//initialize matrix class
	matrix<Fp>::moduleOper = &trunc_oper.FpMod_opers;
	matrix<Gamma>::moduleOper = &trunc_oper.GammaMod_opers;

	//initialize the comod to a trivial one with one generator at degree 0
	trunc_oper.set_to_trivial(comod,0);

	//initialize the curtis tables
	curtis_table<Fp>::ModOper = &trunc_oper.FpMod_opers;
	ResolutionTables.resize(resolution_length+2);
	for(int i=0; i<resolution_length+2; ++i)
		resolutionTables.push_back(&ResolutionTables[i]);
}

//initialize a generic comodule
TruncComodInit::TruncComodInit() : TruncCoMod_generic(&coaction_matrix){
}

//do resolutions
void TruncInit::resolve(string director, std::vector<std::vector<int>> *basis_orders){
	trunc_oper.pre_resolution_tab(comod, director + "maps", director + "gens", resolution_length, resolutionTables, gens, &inj, &indj, &qut, &new_map, director + "BPtables");
}

//save the resolution tables
void TruncInit::saveResolutionTables(string table_data){
	std::fstream tables(table_data, std::ios::out | std::ios::binary);
	for(unsigned i=0; i<ResolutionTables.size(); ++i){
		ResolutionTables[i].save(tables);
	}
}

//save the generators, 4 bytes size, then for each of the generator sets
void TruncInit::save_gens(string gens_data){
	std::fstream genfile(gens_data, std::ios::out | std::ios::binary);
	int32_t sz = gens.size();
	genfile.write((char*)&sz, 4);
	for(int i=0; i<sz; ++i){
		int32_t ss = gens[i].size();
		genfile.write((char*)&ss, 4);
		for(int j=0; j<ss; ++j)
			genfile.write((char*)&gens[i][j], 4);
	}
}
