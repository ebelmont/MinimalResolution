//Steenrod_generic_init.cpp
#include"Steenrod_generic_init.h"

void SteenrodGenericInit::set_comodule(int rank, const std::vector<int> &degree,
                                         std::function<vectors<matrix_index,P>(int)> coaction_rows){
	if((int)degree.size() != rank){
		std::cerr << "SteenrodGenericInit::set_comodule: degree.size() (" << degree.size()
		          << ") does not match rank (" << rank << ")\n" << std::flush;
	}

	//sort each row and check indices are in range, same as BPGenericInit::set_comodule
	std::function<vectors<matrix_index,P>(int)> checked_rows =
	    [rank, &coaction_rows](int i) -> vectors<matrix_index,P> {
		auto row = coaction_rows(i);
		for(auto &tm : row.dataArray){
			if((int)tm.ind < 0 || (int)tm.ind >= rank){
				std::cerr << "SteenrodGenericInit::set_comodule: coaction row " << i
				          << " has an out-of-range index " << tm.ind
				          << " (rank is " << rank << ")\n" << std::flush;
			}
		}
		row.sort();
		return row;
	};

	comod.base_module.rank = rank;
	comod.base_module.degree = degree;
	//ComodInit shadows the inherited (public) coaction_matrix pointer with its
	//own (private) matrix_mem<P> object of the same name -- qualify explicitly
	//to reach comodule_generic's public pointer, same issue as BPComodInit
	//(see BP_generic_init.cpp)
	comod.SteenrodCoMod_generic::coaction_matrix->construct(rank, checked_rows);
}
