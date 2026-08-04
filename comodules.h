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
};

//Every comodule the driver knows about.
const std::vector<ComoduleSpec>& all_comodules();

//The comodule used when none is named on the command line (the sphere).
string default_comodule_name();

//Look a comodule up by name. Returns NULL if there is no such comodule.
const ComoduleSpec* find_comodule(string const &name);

//A printable table of the available comodules, for --list and usage errors.
string list_comodules();
