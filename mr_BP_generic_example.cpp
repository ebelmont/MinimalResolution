//mr_BP_generic_example.cpp
//
//Template + worked example for computing the algebraic Novikov E2 page of a
//user-supplied BP_*BP-comodule (rather than the sphere). To use this for a
//real finite complex:
//
//  1. Copy this file to a new name (e.g. mr_BP_myComplex.cpp) and add it to
//     a copy of BP_generic_compile in place of this file.
//  2. Edit build_comodule() below to describe YOUR complex: its rank (number
//     of BP_*-module generators, i.e. cells/summands), each generator's
//     degree, and its coaction (see the comments inside build_comodule()).
//  3. Rebuild and run exactly like mr_BP: ./mr_BP_myComplex <halfT> <resolution_length>
//     -- it still needs a matching `BPtab <halfT>` run first, for the same
//     reason mr_BP does (the BP_*BP structure maps, not your comodule, come
//     from there). See docs/BUILD_AND_RUN.md.
//
//As shipped, build_comodule() below reconstructs the TRIVIAL comodule (i.e.
//BP_* itself, rank 1, coaction = "1 times itself"). This is deliberate:
//since resolving the trivial comodule is exactly what mr_BP already does,
//running this example and comparing its algNov/Boc output to a normal
//mr_BP run (for the same <halfT>/<resolution_length>) is a correctness
//check on this new code path itself -- if they disagree, something is
//wrong with this new code, not with your eventual real complex.
//
//HOW THIS WORKS: computing a resolution directly over BP_* doesn't work
//(see BP_generic_init.h's comments for why) -- so this program first
//reduces your comodule mod I = (p,v1,v2,...) to get a comodule over the
//FIELD P = BP_*BP/I, resolves THAT with the classical (Steenrod-style)
//machinery, and then lifts that resolution to a genuine BP_*BP-comodule
//resolution of your comodule. You only write your comodule's data once, in
//terms of BP_*BP -- BP_mod_I.h's reduce_coaction_rows_mod_I() derives the
//P-side reduction automatically.
//
//IMPORTANT: nothing here verifies that a hand-supplied coaction actually
//satisfies the comodule axioms (coassociativity, counitality) -- see
//BP_generic_init.h's comments. Get this right independently before trusting
//the output for a real complex.

#include"BP_generic_init.h"
#include"Steenrod_generic_init.h"
#include"BP_mod_I.h"

//EDIT THIS FUNCTION to describe your own complex. BP_oper is the already-
//constructed BP_*BP Hopf algebroid (loaded from a BPtab run) -- use
//BP_oper.BPBP_opers for its ring operations to build BP_*BP elements.
void build_comodule(BP_Op &BP_oper, int &rank, std::vector<int> &degree,
                     std::function<vectors<matrix_index,BPBP>(int)> &coaction_rows){
	//--- EDIT: how many BP_*-module generators does your comodule have? ---
	rank = 1;

	//--- EDIT: the internal degree of each generator ---
	degree = {0};

	//--- EDIT: the coaction on each generator. coaction_rows(i) must return
	//    generator i's coaction as a sparse vectors<matrix_index,BPBP>: a
	//    list of (j, c) pairs, j the index (0..rank-1) of another generator,
	//    c the BP_*BP element multiplying it. Build c with BP_oper's own
	//    ring operations (BP_oper.BPBP_opers is a BPBP_Op -- a RingOp<BPBP>
	//    plus polynomial operations), e.g.:
	//      BP_oper.BPBP_opers.unit(1)         -- the element "1" of BP_*BP
	//      BP_oper.BPBP_opers.monomial(e, c)  -- a single monomial, see WARNING below
	//      BP_oper.BPBP_opers.add(x, y)       -- x + y
	//      BP_oper.BPBP_opers.multiply(x, y)  -- x * y
	//
	//    WARNING -- BPBP's two exponent slots are the opposite of what the
	//    type name suggests. Per BP.cpp:69 ("the right unit, vn is in the
	//    outer"), in BPBP = polynomial<BP> the OUTER exponent indexes the
	//    v_i (via the right unit eta_R), and the INNER coefficient BP's
	//    exponent indexes the t_i. So
	//        BPBP_opers.monomial(singleVar(1,1), unit(1))
	//    is eta_R(v_1), NOT t_1. To build t_1:
	//        BP   inner = BP_oper.monomial(singleVar(1,1), BP_oper.Z3_oper->unit(1));
	//        BPBP t1    = BP_oper.BPBP_opers.monomial(0, inner);
	//    or, simpler and impossible to get backwards, just use BP_oper.h0(),
	//    which computes (eta_R(v1)-eta_L(v1))/p = t_1 from the loaded tables.
	//    Getting this backwards fails SILENTLY: eta_R(v_1) is not in the
	//    augmentation ideal, so the resulting coaction violates counitality,
	//    and nothing here checks that.
	//    As shipped: generator 0's coaction is "1 times itself" -- the
	//    trivial comodule -- exactly what BP_Op::set_to_trivial builds
	//    (compare hopf_algebroid/12.h:1-12).
	coaction_rows = [&BP_oper](int i) -> vectors<matrix_index,BPBP> {
		vectors<matrix_index,BPBP> row;
		BPBP one = BP_oper.BPBP_opers.unit(1);
		row.push({(matrix_index)i, one});
		return row;
	};
}

int main(int argc, char** argv){
	int max_degree = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);

	string filename0 = string(argv[1]) + "_";
	string bp_dir = filename0 + "gBP";       //this program's own BP-side output prefix
	string model_dir = filename0 + "gP";     //this program's own P-side (model) output prefix

	//--- construct the BP-side driver first, so build_comodule() can use its
	//    BP_oper to build BP_*BP elements. Loads the structure maps a BPtab
	//    run already produced for this <halfT> -- same requirement as mr_BP.
	BPGenericInit BPoper(max_degree, resolution_length, filename0+"etaL", filename0+"R2L", filename0+"delta", bp_dir);

	int rank;
	std::vector<int> degree;
	std::function<vectors<matrix_index,BPBP>(int)> coaction_rows;
	build_comodule(BPoper.BP_oper, rank, degree, coaction_rows);

	BPoper.set_comodule(rank, degree, coaction_rows);

	//--- phase 1: resolve the SAME comodule's reduction mod I, as a
	//    P-comodule, with the classical (field-based, correct-as-is)
	//    machinery. The model's resolution_length needs to be at least one
	//    more than the BP-side length, mirroring the mr_st/mr_BP convention
	//    (see docs/BUILD_AND_RUN.md).
	SteenrodGenericInit stOper(3, max_degree, resolution_length+1, model_dir+"steenrod_coaction.data");
	stOper.set_comodule(rank, degree, reduce_coaction_rows_mod_I(coaction_rows));
	stOper.resolve(model_dir);
	stOper.save_gens(model_dir + "gens_data");

	//--- phase 2: lift that model resolution into a genuine BP_*BP-comodule
	//    resolution of your comodule
	std::cout << "starting resolution..." << std::flush;
	BPoper.resolve(model_dir + "gens_data", model_dir + "BPtables");

	//construct the resolution (same post-processing pass mr_BP uses)
	BPoper.resolution();

	//compute algebraic Novikov and Bockstein pages (unchanged from mr_BP --
	//these only ever read back the resolution's own output files)
	BPoper.make_algNov();
	BPoper.make_Boc();

	return 0;
}
