//comodules.h
//
//A registry of BP_*BP-comodules that mr_BP_comod can be asked to resolve,
//selected by name on the command line. The sphere (the trivial comodule
//BP_* itself) is the default, so `./mr_BP_comod <halfT> <len>` reproduces
//what the shipped mr_BP already computes.
//
//TO ADD YOUR OWN COMODULE: see the clearly marked section at the bottom of
//comodules.cpp. It is two steps -- write a builder function, then add one
//row to the table -- and nothing else in the program needs to change.
#pragma once
#include"BP.h"
#include<vector>
#include<functional>

//The signature every comodule builder has. Given the already-initialized
//Hopf algebroid BP_oper (its structure tables are loaded by the time this is
//called, so BP_oper.h0() and friends are usable), fill in:
//  rank          -- the number of BP_*-module generators (e.g. cells)
//  degree        -- degree[i] is generator i's internal degree; size == rank.
//                   These are FULL topological degrees, the same units
//                   exponents.cpp's xnDegs uses: |v_n| = |t_n| = 2(3^n - 1),
//                   so |v_1| = |t_1| = 4 at p=3.
//  coaction_rows -- coaction_rows(i) returns generator i's coaction as a
//                   sparse list of (j, c) pairs: j an index in [0,rank),
//                   c the BP_*BP element multiplying generator j. Rows need
//                   not be sorted; set_comodule sorts them.
//
//The builder's signature is unchanged by the height: the coaction of a
//height-n comodule is written in BP_*BP/I_n, whose elements are BPBP values
//like any other (BP_oper is already reducing mod I_n by the time the builder
//runs, so BP_oper.h0() and friends hand back the reduced elements, and
//anything the builder writes is reduced again by set_comodule).
typedef void (*ComoduleBuilder)(BP_Op &BP_oper, int &rank,
                                 std::vector<int> &degree,
                                 std::function<vectors<matrix_index,BPBP>(int)> &coaction_rows);

//One entry in the registry.
struct ComoduleSpec{
	//the name typed on the command line
	string name;
	//one-line summary, shown by --list and in the usage message
	string description;
	//builds the comodule's rank/degrees/coaction
	ComoduleBuilder build;
	//The height n: the comodule's underlying module is free over
	//BP_*/I_n, where I_n = (p, v_1, ..., v_{n-1}).
	//
	//  0 -- free over BP_* itself. The classical case.
	//  1 -- free over BP_*/p        (e.g. BP_*(S/p))
	//  2 -- free over BP_*/(p,v_1)  (e.g. BP_*(S/(p,v_1)) = BP_*(V(1)))
	//
	//I_n is an invariant ideal, so (BP_*/I_n, BP_*BP/I_n) is again a Hopf
	//algebroid, and for a comodule M with I_n M = 0 the change-of-rings
	//isomorphism
	//    Ext_{BP_*BP}(BP_*, M) = Ext_{BP_*BP/I_n}(BP_*/I_n, M)
	//says resolving M over the quotient computes the same Ext. The whole
	//pipeline then runs over BP_*/I_n: see BP_Op::height.
	//
	//n is a property of the comodule, not a free parameter, so it is
	//recorded here rather than asked for on the command line.
	//
	//GETTING IT WRONG IS SILENT IN BOTH DIRECTIONS, and nothing checks it:
	//  too small -- the machinery treats a module that is not free over
	//               BP_*/I_n as though it were, and the answer is nonsense;
	//  too large -- the run is internally consistent but computes
	//               Ext(BP_*, M/I_n) instead of Ext(BP_*, M). Declaring the
	//               sphere at n=1, for instance, quietly computes the E2
	//               page of S/p.
	//State n from a proof that the module is free over BP_*/I_n and over
	//nothing larger, not from a guess.
	int height;
};

//Every comodule the driver knows about.
const std::vector<ComoduleSpec>& all_comodules();

//The comodule used when none is named on the command line (the sphere).
string default_comodule_name();

//Look a comodule up by name. Returns NULL if there is no such comodule.
const ComoduleSpec* find_comodule(string const &name);

//A printable table of the available comodules, for --list and usage errors.
string list_comodules();
