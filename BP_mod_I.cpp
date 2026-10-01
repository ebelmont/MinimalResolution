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
		//BPBP's OUTER exponent indexes the v_i and its INNER (coefficient)
		//BP's exponent indexes the t_i -- the opposite of what the type name
		//suggests; see the warning at the top of comodules.cpp. So a term
		//with a nonzero outer exponent carries a v, hence lies in
		//I = (p,v_1,v_2,...) and dies. (This holds whichever unit the outer
		//slot uses: eta_L(v_n) is in I by definition, and eta_R(v_n) is in
		//I*BPBP because I is invariant.)
		if(tm.ind != 0) continue;
		//what is left is a polynomial in the t_i with Z_3 coefficients;
		//reduce those mod p. The inner terms are already in index order, so
		//pushing them in order keeps the result canonical.
		for(auto &im : tm.coeficient.dataArray){
			Fp c = (Fp)(im.coeficient % 3);
			if(c!=0)
				result.push({im.ind, c});
		}
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
