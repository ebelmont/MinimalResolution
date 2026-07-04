//BP_mod_I.h
#pragma once
#include"BP.h"
#include"steenrod.h"

//Reduce an element of BP (= BP_*, a polynomial in v_1,v_2,... over Z3) to
//its class in BP_*/(p, v_1, v_2, ...) = F_3, i.e. its constant term (v_i-
//exponent 0), reduced mod 3. Every term with a nonzero v_i-exponent maps to
//0, since it lies in the ideal I = (p, v_1, v_2, ...).
Fp reduce_BP_mod_I(BP const &x);

//Reduce an element of BPBP (= BP_*BP, a polynomial in t_1,t_2,... with BP
//coefficients) to its class in P = BP_*BP/I (the same polynomial in t_i,
//each BP coefficient reduced mod I via reduce_BP_mod_I above). This relies
//on BPBP and P sharing the same monomial (t_i) encoding -- true here because
//both are built from the same shared exponents.h/exponents.cpp table (see
//docs/pipelines/STEENROD.md sec 1 and docs/ARCHITECTURE.md sec 1 for why
//BP_*BP/I's t_i and the classical pipeline's generators share a grading and
//encoding by construction, not coincidence).
P reduce_BPBP_mod_I(BPBP const &x);

//Reduce one coaction-matrix row (a sparse vectors<matrix_index,BPBP>) to its
//class mod I (a sparse vectors<matrix_index,P>), dropping any entry whose
//reduction is zero.
vectors<matrix_index,P> reduce_row_mod_I(vectors<matrix_index,BPBP> const &row);

//Wrap a BP_*BP-valued coaction-row function (as you would pass to
//BPGenericInit::set_comodule) into the corresponding P-valued one (as you
//would pass to SteenrodGenericInit::set_comodule), by reducing every row mod
//I. Since M is free over BP_*, M/I has the SAME rank and generator degrees
//as M -- only the coaction needs reducing. This is the one function that
//lets you enter a comodule's data once (in BP_*BP) and get both phases of
//the computation (see mr_BP_generic_example.cpp) from it.
std::function<vectors<matrix_index,P>(int)> reduce_coaction_rows_mod_I(
    std::function<vectors<matrix_index,BPBP>(int)> bp_rows);
