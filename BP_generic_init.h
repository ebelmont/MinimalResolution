//BP_generic_init.h
#pragma once
#include"BP_init.h"

//BPGenericInit resolves a *user-supplied* BP_*BP-comodule (rank, per-
//generator degree, and coaction), instead of the trivial rank-1 comodule
//(i.e. BP_* itself) that BPInit's constructor sets up via
//BP_Op::set_to_trivial. Use this to compute the algebraic Novikov E2 page of
//any finitely generated BP_*BP-comodule that is free over BP_* -- e.g. the
//BP-homology of a finite complex -- rather than of the sphere.
//
//This still needs the BP Hopf algebroid's own structure maps (etaL/delta/
//R2L), exactly as the sphere computation does -- those describe (BP_*,
//BP_*BP) itself, not the comodule being resolved, so they're unchanged and
//still come from a BPtab run.
//
//IMPORTANT -- read this before using resolve(): there is no way to search
//for a resolution of a BP_*BP-comodule directly (Gaussian-elimination-style
//pivot search, as Hopf_Algebroid::pre_resolution_tab does, relies on
//"invertible leading term," and PolyOp::invertible only recognizes bare
//constants in BP_* = Z3[v1,v2,...] as invertible -- almost nothing a real
//coaction produces qualifies). This is why the sphere computation
//(BPInit::resolve()) never calls pre_resolution_tab with ring=BP either --
//it exclusively lifts a resolution computed over the FIELD P = BP_*BP/I
//instead. resolve() below does the same: it takes an already-computed model
//resolution of this SAME comodule's reduction mod I (built with
//SteenrodGenericInit, see Steenrod_generic_init.h and BP_mod_I.h) and lifts
//it via Hopf_Algebroid::pre_resolution_modeled. See
//comodules.cpp / mr_BP_comod.cpp for how the two phases fit together, and
//docs/CODE_WALKTHROUGH.md / docs/ARCHITECTURE.md section 5 for the
//background this assumes.
//
//ALSO IMPORTANT: set_comodule() does not, and cannot easily, verify that the
//coaction you supply actually satisfies the comodule axioms (coassociativity
//with the algebroid's comultiplication, and counitality). A mistake here
//will not crash -- it will silently produce a wrong Ext computation. It is
//your responsibility to make sure the coaction matrix you supply is
//mathematically correct; this class only checks that the data you give it
//is structurally well-formed (indices in range, rows sorted).
class BPGenericInit : public BPInit{
public:
	//re-use BPInit's constructor: this still needs max_deg, resolution_length,
	//and the BPtab-produced etaL/delta/R2L structure-map files, exactly as
	//for the sphere computation
	using BPInit::BPInit;

	//Populate the comodule to be resolved, replacing the trivial rank-1
	//comodule BPInit's constructor set up.
	//  rank    -- the number of BP_*-module generators of the comodule
	//  degree  -- degree[i] is generator i's internal grading (size must be rank)
	//  coaction_rows -- coaction_rows(i) must return generator i's coaction,
	//    as a sparse vectors<matrix_index,BPBP>: a list of (j, c) pairs where
	//    j ranges over the OTHER generators 0..rank-1 and c is a BP_*BP
	//    element (build these with BP_oper.BPBP_opers, e.g. .unit(1) for the
	//    element "1", .monomial(...)/.constant(...) otherwise). Each
	//    returned vectors<matrix_index,BPBP> need not already be sorted by
	//    index -- set_comodule sorts it for you -- but every index it uses
	//    must be in range [0, rank).
	//See comodules.cpp for complete worked examples.
	void set_comodule(int rank, const std::vector<int> &degree,
	                   std::function<vectors<matrix_index,BPBP>(int)> coaction_rows);

	//Lift the comodule set by set_comodule's resolution from an
	//already-computed model resolution of its reduction mod I (built with
	//SteenrodGenericInit -- see BP_mod_I.h's reduce_coaction_rows_mod_I()
	//for deriving that reduction automatically from the same BP_*BP-valued
	//coaction passed to set_comodule).
	//  model_gens_data  -- the file written by the model's
	//                      SteenrodGenericInit::save_gens()
	//  model_table_file -- the file written by the model's
	//                      SteenrodGenericInit::resolve()'s tablename argument
	//                      (i.e. "<director>BPtables" for whatever director
	//                      you resolved the model with)
	//Writes <director>maps / <director>gens exactly as BPInit::resolve()
	//does, so BPInit::resolution() / make_algNov() / make_Boc() /
	//mult_table() all work unchanged afterward -- they only ever read those
	//files back, never anything about which comodule produced them.
	void resolve(string model_gens_data, string model_table_file);
};
