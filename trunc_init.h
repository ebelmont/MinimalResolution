//trunc_init.h
#pragma once
#include"trunc_hopf.h"
#include"matrices_mem.h"

//initialize a generic comodule
class TruncComodInit : public TruncCoMod_generic{
	matrix_mem<Gamma> coaction_matrix;
public:
	//the constructor
	TruncComodInit();
};

//the initialization of the data for Gamma = F_p[x]/(x^GammaDim - 1)
class TruncInit{
public:
	//the length of the resolution
	int resolution_length;

	//the operations on the Hopf algebroid (F_p, Gamma)
	TruncHopf_Op trunc_oper;

	//the matrix for the co-multiplication
	matrix_mem<Gamma> deltaTable;

	//the container for the matrix of the injection to a cofree comodule, and the quotient to the next comodule
	matrix_mem<Fp> inj, indj, qut, new_map;

	//the curtis table for the resolutions
	std::vector<curtis_table<Fp>*> resolutionTables;
	std::vector<curtis_table_mem<Fp>> ResolutionTables;

	//the generators for a resolution
	std::vector<std::vector<int>> gens;

	//the comodule, initialized to the trivial comodule of rank 1 (F_p as trivial comodule)
	TruncComodInit comod;

	//the constructor
	TruncInit(int prime, int max_deg, int resolution_length);

	//do resolutions
	void resolve(string director, std::vector<std::vector<int>> *basis_orders = NULL);

	//save the resolution tables
	void saveResolutionTables(string);

	//save the generators
	void save_gens(string);
};
