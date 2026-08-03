//mr_BP_S_alpha1.cpp
//
//Computes the algebraic Novikov E2 page for BP_*(S/alpha_1) at p=3, where
//alpha_1 in pi_3(S)_(3) = pi_{2p-3}(S)_(p) is the first element of the alpha
//family and S/alpha_1 = cofib(S^3 --alpha_1--> S^0) = S^0 ∪_{alpha_1} e^4.
//
//THE COMODULE (see docs/GENERAL_COMODULES.md for how to run this):
//
//  BP_* is concentrated in even degrees, so BP_3 = 0 and hence (alpha_1)_* = 0
//  on BP-homology. The cofiber long exact sequence therefore collapses to a
//  short exact sequence of BP_*BP-comodules
//
//      0 --> BP_* --> BP_*(S/alpha_1) --> Sigma^4 BP_* --> 0
//
//  so BP_*(S/alpha_1) is FREE over BP_* of rank 2, on a bottom-cell class x_0
//  in degree 0 and a top-cell class x_4 in degree 4 (|t_1| = |v_1| = 2(p-1) = 4
//  at p=3 -- the same units exponents.cpp's xnDegs uses).
//
//  The bottom cell is a sub-comodule, so its coaction is trivial. The extension
//  class of the sequence above lives in Ext^1_{BP_*BP}(Sigma^4 BP_*, BP_*) =
//  Ext^{1,4}_{BP_*BP}(BP_*,BP_*) = Z/3, and it is by construction the class
//  detecting alpha_1 -- namely [t_1], the generator usually written h_0.
//  So the coaction matrix is
//
//      psi(x_0) = 1 (x) x_0
//      psi(x_4) = 1 (x) x_4  +  t_1 (x) x_0
//
//  i.e., writing the coaction as a matrix acting on (x_0, x_4):
//
//         [  1    0  ]
//         [ t_1   1  ]
//
//  Coassociativity holds because t_1 is primitive (Delta(t_1) = t_1(x)1 + 1(x)t_1),
//  and counitality because eps(t_1) = 0.
//
//  Rather than hand-building t_1, this file uses BP_Op::h0(), which the repo
//  already provides: it computes (eta_R(v_1) - eta_L(v_1))/p from the loaded
//  structure tables, which at p=3 is exactly t_1 (eta_R(v_1) = v_1 + p*t_1 for
//  Hazewinkel generators). Using it means the element is guaranteed consistent
//  with this codebase's own conventions, and h0() prints itself so you can
//  eyeball that it really is t_1 before trusting the run.
//
//NOTE ON WHICH alpha_1: this is the cofiber of the element alpha_1 of pi_3(S),
//a map between SPHERES. Do not confuse it with the Adams self map
//Sigma^{2p-2} V(0) -> V(0) (whose cofiber is the Smith-Toda complex V(1) =
//S/(p,v_1)); that is a different complex with a different, rank-4 comodule.
//
//Run exactly like mr_BP, after a matching `BPtab <halfT>` run:
//    ./mr_BP_S_alpha1 <halfT> <resolution_length>
//<halfT> must be at least 4 for t_1 to exist at all; use something comfortably
//larger (e.g. 30+) for a resolution that shows anything interesting.

#include"BP_generic_init.h"
#include"Steenrod_generic_init.h"
#include"BP_mod_I.h"

void build_comodule(BP_Op &BP_oper, int &rank, std::vector<int> &degree,
                     std::function<vectors<matrix_index,BPBP>(int)> &coaction_rows){
	//rank 2: the bottom cell (degree 0) and the top cell (degree 4)
	rank = 2;
	degree = {0, 4};

	//t_1, computed from the structure tables as (eta_R(v_1)-eta_L(v_1))/p.
	//Captured by value so it is computed once, not once per coaction row.
	BPBP t1 = BP_oper.h0();

	coaction_rows = [&BP_oper, t1](int i) -> vectors<matrix_index,BPBP> {
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

int main(int argc, char** argv){
	int max_degree = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);

	string filename0 = string(argv[1]) + "_";
	string bp_dir = filename0 + "a1BP";      //this program's own BP-side output prefix
	string model_dir = filename0 + "a1P";    //this program's own P-side (model) output prefix

	BPGenericInit BPoper(max_degree, resolution_length, filename0+"etaL", filename0+"R2L", filename0+"delta", bp_dir);

	int rank;
	std::vector<int> degree;
	std::function<vectors<matrix_index,BPBP>(int)> coaction_rows;
	build_comodule(BPoper.BP_oper, rank, degree, coaction_rows);

	BPoper.set_comodule(rank, degree, coaction_rows);

	//phase 1: resolve the reduction mod I = (p,v_1,v_2,...) over the field-based
	//Hopf algebra P = BP_*BP/I. Note t_1 does NOT die mod I, so the reduced
	//comodule is still the nonsplit extension -- this is the "algebraic alpha_1".
	SteenrodGenericInit stOper(3, max_degree, resolution_length+1, model_dir+"steenrod_coaction.data");
	stOper.set_comodule(rank, degree, reduce_coaction_rows_mod_I(coaction_rows));
	stOper.resolve(model_dir);
	stOper.save_gens(model_dir + "gens_data");

	//phase 2: lift that model resolution to BP_*BP
	std::cout << "starting resolution..." << std::flush;
	BPoper.resolve(model_dir + "gens_data", model_dir + "BPtables");

	BPoper.resolution();
	BPoper.make_algNov();
	BPoper.make_Boc();

	return 0;
}
