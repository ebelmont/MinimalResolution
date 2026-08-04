//BP_generic_init.cpp
#include"BP_generic_init.h"

void BPGenericInit::set_comodule(int rank, const std::vector<int> &degree,
                                   std::function<vectors<matrix_index,BPBP>(int)> coaction_rows){
	if((int)degree.size() != rank){
		std::cerr << "BPGenericInit::set_comodule: degree.size() (" << degree.size()
		          << ") does not match rank (" << rank << ")\n" << std::flush;
	}

	//wrap the caller's row function so every row is sorted by index (as
	//every sparse vectors<matrix_index,R> elsewhere in this codebase is
	//assumed to be) and so out-of-range indices are caught early, before
	//they cause a confusing failure deep inside the resolution algorithm
	std::function<vectors<matrix_index,BPBP>(int)> checked_rows =
	    [rank, &coaction_rows](int i) -> vectors<matrix_index,BPBP> {
		auto row = coaction_rows(i);
		for(auto &tm : row.dataArray){
			if((int)tm.ind < 0 || (int)tm.ind >= rank){
				std::cerr << "BPGenericInit::set_comodule: coaction row " << i
				          << " has an out-of-range index " << tm.ind
				          << " (rank is " << rank << ")\n" << std::flush;
			}
		}
		row.sort();
		return row;
	};

	comod.base_module.rank = rank;
	comod.base_module.degree = degree;
	//BPComodInit shadows the inherited (public) coaction_matrix pointer with
	//its own (private) matrix_file<BPBP> object of the same name -- qualify
	//explicitly to reach comodule_generic's public pointer, which already
	//points at that same object (see BPComodInit's constructor, BP_init.cpp:5-6)
	comod.BPCoMod_generic::coaction_matrix->construct(rank, checked_rows);
}

void BPGenericInit::resolve(string model_gens_data, string model_table_file){
	//load the model's chosen generator indices per step (same file format
	//BPInit::load_gens already reads for the sphere computation, just
	//pointed at our own model instead of mr_st's)
	load_gens(model_gens_data);

	//lift F_p-valued vectors from the model into BP-valued ones -- the same
	//lift BPInit::resolve() uses for the sphere, unchanged: it's a property
	//of the Hopf algebroid (BP_*, BP_*BP), not of which comodule is resolved
	std::function<vectors<matrix_index,BP>(const vectors<matrix_index,Fp>&)> tfm = [this] (const vectors<matrix_index,Fp>& v){
		return BP_oper.lift(v); };

	BP_oper.pre_resolution_modeled(comod, director+"maps", director+"gens",
	                                 resolution_length, model_table_file, &ctable, gens, tfm,
	                                 &inj, &qut, &indj, &new_map, director+"back");

	//combine the per-step generator files into one, exactly as BPInit::resolve() does
	BP_oper.gens_file_combiner(director + "gens", resolution_length, comod);
}
