//BP_mod_I.cpp
#include"BP_mod_I.h"

Fp reduce_BP_mod_I(BP const &x){
	if(x.size()>0 && x.dataArray[0].ind==0)
		return (Fp)(x.dataArray[0].coeficient % 3);
	return 0;
}

P reduce_BPBP_mod_I(BPBP const &x){
	P result;
	for(auto &tm : x.dataArray){
		Fp c = reduce_BP_mod_I(tm.coeficient);
		if(c!=0)
			result.push({tm.ind, c});
	}
	return result;
}

vectors<matrix_index,P> reduce_row_mod_I(vectors<matrix_index,BPBP> const &row){
	vectors<matrix_index,P> result;
	for(auto &tm : row.dataArray){
		P c = reduce_BPBP_mod_I(tm.coeficient);
		if(c.size()>0)
			result.push({tm.ind, c});
	}
	return result;
}

std::function<vectors<matrix_index,P>(int)> reduce_coaction_rows_mod_I(
    std::function<vectors<matrix_index,BPBP>(int)> bp_rows){
	return [bp_rows](int i) -> vectors<matrix_index,P>{
		return reduce_row_mod_I(bp_rows(i));
	};
}
