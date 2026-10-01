//BP_mod_I.h
#pragma once
#include"BP.h"
#include"steenrod.h"

//Reduce an element of BP (= BP_*, a polynomial in v_1,v_2,... over Z3) to
//its class in BP_*/(p, v_1, v_2, ...) = F_3, i.e. its constant term (v_i-
//exponent 0), reduced mod 3. Every term with a nonzero v_i-exponent maps to
//0, since it lies in the ideal I = (p, v_1, v_2, ...).
//
//NOTE this is for genuine BP_* elements only. It is deliberately NOT used by
//reduce_BPBP_mod_I below: the BP sitting in BPBP's coefficient slot is a
//polynomial in the t_i, not the v_i, so "kill every term with a nonzero
//exponent" is the wrong operation there.
Fp reduce_BP_mod_I(BP const &x);

//Reduce an element of BPBP = BP_*BP to its class in P = BP_*BP/I.
//
//MIND THE SLOTS. BPBP = polynomial<BP> is NOT "a polynomial in t_1,t_2,...
//with BP_* coefficients" in the layout the type name suggests: the OUTER
//exponent indexes the v_i and the INNER (coefficient) BP's exponent indexes
//the t_i -- see the warning at the top of comodules.cpp, and BP.cpp:69. So
//the reduction keeps only the terms whose OUTER exponent is 0 (any v makes
//the term lie in I, whichever unit put it there: eta_L(v_n) is in I by
//definition and eta_R(v_n) is in I*BPBP because I is invariant), and reduces
//the remaining inner t-polynomial's Z_3 coefficients mod p.
//
//Two values pin the convention down, and are worth re-checking after any
//change here: t_1 (= BP_Op::h0()) must reduce to t_1, and eta_R(v_1)
//(= BPBP_opers.monomial(singleVar(1,1), unit(1))) must reduce to 0. Getting
//these backwards is silent -- it reduces every off-diagonal coaction entry
//to 0, i.e. it hands the second phase a SPLIT comodule.
//
//The result relies on BPBP and P sharing the same monomial (t_i) encoding --
//true here because both are built from the same shared exponents.h/
//exponents.cpp table (see docs/pipelines/STEENROD.md sec 1 and
//docs/ARCHITECTURE.md sec 1 for why BP_*BP/I's t_i and the classical
//pipeline's generators share a grading and encoding by construction, not
//coincidence).
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
//the computation (see mr_BP_comod.cpp) from it.
std::function<vectors<matrix_index,P>(int)> reduce_coaction_rows_mod_I(
    std::function<vectors<matrix_index,BPBP>(int)> bp_rows);
