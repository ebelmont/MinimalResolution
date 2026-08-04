//comodules.cpp
//
//The comodules mr_BP_comod can resolve, and the table that maps command-line
//names to them. See the "ADDING YOUR OWN COMODULE" section at the bottom.
//
//=====================================================================
//WARNING -- BPBP's two exponent slots are the opposite of what the type
//name suggests. Read this before writing any coaction by hand.
//
//BPBP = polynomial<BP> is NOT "polynomial in the t_i with BP_* coefficients"
//in the layout you would guess. Per BP.cpp:69 ("the right unit, vn is in the
//outer"):
//    * the OUTER exponent indexes the v_i, included via the RIGHT unit eta_R
//    * the INNER (coefficient) BP's exponent indexes the t_i
//So the obvious-looking
//        BP_oper.BPBP_opers.monomial(singleVar(1,1), BP_oper.unit(1))
//is eta_R(v_1), NOT t_1. The correct t_1 is
//        BP   inner = BP_oper.monomial(singleVar(1,1), BP_oper.Z3_oper->unit(1));
//        BPBP t1    = BP_oper.BPBP_opers.monomial(0, inner);
//or, better, just call BP_oper.h0(), which computes
//(eta_R(v1) - eta_L(v1))/p = t_1 from the loaded structure tables and cannot
//be gotten backwards. (Verified: h0() is bitwise equal to the construction
//above, and both differ from the monomial(singleVar(1,1),...) form.)
//
//Getting this backwards fails SILENTLY: eta_R(v_1) is not in the
//augmentation ideal, so the resulting coaction violates counitality -- and
//nothing in this pipeline checks the comodule axioms.
//=====================================================================
#include"comodules.h"

//---------------------------------------------------------------------
//the sphere: the trivial comodule BP_* itself
//---------------------------------------------------------------------
//Rank 1, one generator in degree 0, coaction "1 times itself" -- exactly
//what BP_Op::set_to_trivial builds (hopf_algebroid/12.h:1-12), which is what
//the shipped mr_BP resolves. Running this comodule through the generic path
//and comparing its algNov/Boc tables against a plain mr_st+BPtab+mr_BP run
//is therefore a real regression test of the generic machinery itself; see
//docs/GENERAL_COMODULES.md.
static void build_sphere(BP_Op &BP_oper, int &rank, std::vector<int> &degree,
                          std::function<vectors<matrix_index,BPBP>(int)> &coaction_rows){
	rank = 1;
	degree = {0};

	coaction_rows = [&BP_oper](int i) -> vectors<matrix_index,BPBP>{
		vectors<matrix_index,BPBP> row;
		row.push({(matrix_index)i, BP_oper.BPBP_opers.unit(1)});
		return row;
	};
}

//---------------------------------------------------------------------
//S/alpha_1: the cofiber of alpha_1 in pi_3(S)_(3)
//---------------------------------------------------------------------
//alpha_1 in pi_{2p-3}(S)_(p) = pi_3(S)_(3) is the first element of the alpha
//family, and S/alpha_1 = cofib(S^3 --alpha_1--> S^0) = S^0 u_{alpha_1} e^4.
//
//BP_* is concentrated in even degrees, so BP_3 = 0 and hence (alpha_1)_* = 0
//on BP-homology. The cofiber long exact sequence therefore collapses to a
//short exact sequence of BP_*BP-comodules
//
//    0 --> BP_* --> BP_*(S/alpha_1) --> Sigma^4 BP_* --> 0,
//
//so BP_*(S/alpha_1) is FREE over BP_* of rank 2: a bottom-cell class x_0 in
//degree 0 and a top-cell class x_4 in degree 4 (|t_1| = |v_1| = 2(p-1) = 4).
//
//The bottom cell is a sub-comodule, hence primitive. The extension class of
//the sequence lives in Ext^1_{BP_*BP}(Sigma^4 BP_*, BP_*) =
//Ext^{1,4}_{BP_*BP}(BP_*,BP_*) = Z/3, and by the geometric boundary theorem
//it is the class detecting alpha_1 -- namely [t_1], usually written h_0.
//That class is the off-diagonal entry:
//
//    psi(x_0) = 1 (x) x_0
//    psi(x_4) = 1 (x) x_4  +  t_1 (x) x_0
//
//i.e. as a matrix on (x_0, x_4):   [  1    0 ]
//                                  [ t_1   1 ]
//
//Coassociativity holds because t_1 is primitive (Delta t_1 = t_1(x)1 + 1(x)t_1);
//counitality because eps(t_1) = 0. Note the answer is forced up to a unit:
//in degree 4 the reduced part of BP_*BP over BP_* is spanned by t_1 alone,
//so the only real content is that the extension is nonsplit -- which is
//exactly the statement alpha_1 != 0. (Replacing t_1 by 2t_1 rescales x_4 and
//gives an isomorphic comodule.)
//
//NOTE ON WHICH alpha_1: this is the cofiber of an element of pi_3(S), a map
//between SPHERES. It is NOT the cofiber of the Adams self map
//Sigma^4 V(0) -> V(0), whose cofiber is the Smith-Toda complex
//V(1) = S/(3,v_1) -- a different complex with a rank-4 comodule.
static void build_S_alpha1(BP_Op &BP_oper, int &rank, std::vector<int> &degree,
                            std::function<vectors<matrix_index,BPBP>(int)> &coaction_rows){
	rank = 2;
	degree = {0, 4};

	//t_1, taken from the structure tables as (eta_R(v_1)-eta_L(v_1))/p rather
	//than hand-built -- see the WARNING at the top of this file. Computed once
	//here and captured by value, not recomputed per row.
	BPBP t1 = BP_oper.h0();

	coaction_rows = [&BP_oper, t1](int i) -> vectors<matrix_index,BPBP>{
		vectors<matrix_index,BPBP> row;
		BPBP one = BP_oper.BPBP_opers.unit(1);
		if(i == 0){
			//psi(x_0) = 1 (x) x_0 -- the bottom cell is a sub-comodule
			row.push({(matrix_index)0, one});
		} else {
			//psi(x_4) = t_1 (x) x_0 + 1 (x) x_4
			row.push({(matrix_index)0, t1});
			row.push({(matrix_index)1, one});
		}
		return row;
	};
}

//=====================================================================
//ADDING YOUR OWN COMODULE
//
//1. Write a builder function with the same shape as the two above:
//
//      static void build_myComplex(BP_Op &BP_oper, int &rank,
//                                   std::vector<int> &degree,
//                                   std::function<vectors<matrix_index,BPBP>(int)> &coaction_rows){
//          rank   = 2;              //number of BP_*-module generators
//          degree = {0, 4};         //full topological degrees, |v_1|=|t_1|=4
//          coaction_rows = [&BP_oper](int i) -> vectors<matrix_index,BPBP>{
//              vectors<matrix_index,BPBP> row;
//              //...push {j, c} pairs; see build_S_alpha1 above...
//              return row;
//          };
//      }
//
//   Build the BP_*BP elements c with BP_oper.BPBP_opers (.unit(1), .add,
//   .multiply, ...), heeding the WARNING at the top of this file about which
//   exponent slot is which. BP_oper's structure tables are already loaded
//   when your builder runs, so BP_oper.h0() (= t_1) and BP_oper.thetas()
//   (= beta_1 etc.) are available.
//
//2. Add one row to the table below.
//
//That is all -- mr_BP_comod picks the new name up automatically, including
//in --list and in its usage message.
//
//NOTHING HERE CHECKS THE COMODULE AXIOMS. An incorrect coaction will not
//crash; it will silently produce a wrong Ext computation. Verify
//coassociativity and counitality by hand before trusting a run.
//=====================================================================
static const std::vector<ComoduleSpec> comodule_table = {
	{"sphere",    "the sphere: the trivial comodule BP_* (rank 1, degree 0) -- what mr_BP resolves", build_sphere},
	{"alpha_1",   "S/alpha_1 = cofib(S^3 -> S^0) (rank 2, degrees 0 and 4, off-diagonal t_1)",       build_S_alpha1},
};

const std::vector<ComoduleSpec>& all_comodules(){
	return comodule_table; }

string default_comodule_name(){
	return "sphere"; }

const ComoduleSpec* find_comodule(string const &name){
	for(auto const &c : comodule_table)
		if(c.name == name) return &c;
	return NULL;
}

string list_comodules(){
	//pad the names into a column so the descriptions line up
	size_t width = 0;
	for(auto const &c : comodule_table)
		if(c.name.size() > width) width = c.name.size();

	string result;
	for(auto const &c : comodule_table){
		result += "  " + c.name + string(width - c.name.size(), ' ') + "  " + c.description;
		if(c.name == default_comodule_name()) result += "   [default]";
		result += "\n";
	}
	return result;
}
