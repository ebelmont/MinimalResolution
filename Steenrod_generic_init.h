//Steenrod_generic_init.h
#pragma once
#include"steenrod_init.h"

//SteenrodGenericInit resolves a *user-supplied* P-comodule (P = BP_*BP/I,
//rank/degree/coaction), instead of the trivial rank-1 comodule (the sphere's
//own F_3) that SteenrodInit's constructor sets up via set_to_trivial.
//
//This is "phase 1" of computing the algebraic Novikov E2 page of a general
//BP_*BP-comodule M: resolve M/I (a P-comodule, over the FIELD F_3) with this
//class -- unlike the BP side, this is a from-scratch resolution
//(pre_resolution_tab, inherited unchanged from SteenrodInit::resolve()) and
//it is correct as-is, because P's base ring is a field (see BP_generic_init.h
//for why the analogous from-scratch approach does NOT work directly over
//BP_* itself). The resulting resolution is then used as the "model" that
//BPGenericInit lifts into a genuine BP_*BP-comodule resolution of M -- see
//BP_mod_I.h for the M -> M/I reduction that connects the two, and
//mr_BP_comod.cpp for how the two phases are wired together.
class SteenrodGenericInit : public SteenrodInit{
public:
	//re-use SteenrodInit's constructor unchanged: this still needs the
	//prime, max degree, resolution length, and the coproduct-table file
	//(self-generated on first run, exactly as for mr_st)
	using SteenrodInit::SteenrodInit;

	//Populate the P-comodule to be resolved, replacing the trivial rank-1
	//comodule the constructor set up.
	//  rank   -- the number of generators (same rank as the BP_*BP-comodule
	//            M this is the mod-I reduction of, since M is free over BP_*)
	//  degree -- degree[i] is generator i's internal grading (size must be rank)
	//  coaction_rows -- coaction_rows(i) must return generator i's coaction
	//    as a sparse vectors<matrix_index,P>: a list of (j, c) pairs, j the
	//    index (0..rank-1) of another generator, c an element of P. In
	//    practice you will rarely write this by hand -- see
	//    BP_mod_I.h's reduce_coaction_rows_mod_I(), which derives it
	//    automatically from the same BP_*BP-valued coaction you give
	//    BPGenericInit::set_comodule().
	void set_comodule(int rank, const std::vector<int> &degree,
	                   std::function<vectors<matrix_index,P>(int)> coaction_rows);
};
