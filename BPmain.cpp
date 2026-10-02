//BPmain.cpp
#include"BP.h"
#include"BP_init.h"
#include<cstring>
#include<memory>
#include <omp.h>
#include"diag.h"

int main(int argc, char** argv){
	//set the maximal degree
	omp_set_max_active_levels(1);
	int max_degree = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);
	DiagTimer _total("mr_BP");
	//initialize the monomial index
	auto mon_timer = std::unique_ptr<DiagTimer>(new DiagTimer("monomial_index"));
	monomial_index mon_index(max_degree);
	mon_timer.reset();
	
	string filename = argv[1];
	string filename0 = filename + "_";
	filename += "_BP";
	
	//construct the operators for BP
	auto init_timer = std::unique_ptr<DiagTimer>(new DiagTimer("init"));
	BPInit BPoper(max_degree, resolution_length, filename0 + "etaL", filename0 + "R2L", filename0 + "delta", filename);
	
	//load the resolution table for BP/I
//	BPoper.loadResolutionTables(filename0 + "ResTables");
	//load the generators
	BPoper.load_gens(filename0 + "gens_data");
	init_timer.reset();
	
	std::cout << "starting resolution..." << std::flush;
	//do the resolution
	{ DiagTimer _dt("resolve"); BPoper.resolve(); }
	
	//construct the resolution
	{ DiagTimer _dt("resolution"); BPoper.resolution(); }
	
	//compute algebraic Novikov
	{ DiagTimer _dt("make_algNov"); BPoper.make_algNov(); }
	
	//compute Bocstein table
	{ DiagTimer _dt("make_Boc"); BPoper.make_Boc(); }
	
	//compute the multiplicative structure
	{ DiagTimer _dt("mult_table"); BPoper.mult_table(); }
	{ DiagTimer _dt("mult_theta"); BPoper.mult_theta(resolution_length); }
	
	return 0;
}
