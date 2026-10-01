//BP_init.h
#pragma once
#include"BP.h"
#include<memory>
#include"matrices_mem.h"
#include"matrices_stream.h"
#include"algNov.h"
#include"multiplication.h"
#include"Boc.h"

//initialize a generic comodule
class BPComodInit : public BPCoMod_generic{
	matrix_file<BPBP> coaction_matrix;
public:
	//the constructor
	BPComodInit(string);
};

//the initialization of the data for BP
class BPInit{
public:
	//the maximal degee
	int max_degree;
	
	//the length of the resolution
	int resolution_length;
	
	//the director for the data
	string director;
	
	//The height n of the invariant ideal I_n = (p, v_1, ..., v_{n-1}) the
	//whole computation runs modulo; 0 is the classical case (base ring
	//BP_*), which is what mr_BP uses and what every pre-existing caller
	//gets. See BP_Op::height for what n>0 means.
	int height;
	
	//The operations on Z3. Owned through a pointer because the base ring's
	//arithmetic depends on the height: characteristic 0 (Z3_Op) at height 0,
	//characteristic p (Z3_mod_p_Op) above it. Z3_oper is the reference
	//everything else uses, exactly as the plain member did before.
	std::unique_ptr<Z3_Op> Z3_oper_store;
	Z3_Op &Z3_oper;
	//the operations on BP
	BP_Op BP_oper;
	//operations on F3-modules
	ModuleOp<matrix_index,F3> F3Mod_opers;
	
	//the matrices for the structure data
	matrix_file<BP> etaL_matrix, R2L_matrix;
	matrix_file<BPBP> delta_matrix;
	
	//the complex of primitives
	BPComplex Complex;
	//algebraic Novikov tables
	std::vector<algNov_table> AAN_table;
	algNov_tables AANtables;
	//Bockstein tables
	std::vector<Boc_table> B_table;
	Boc_tables Btables;
	
	//the container for the matrix of the injection to a cofree comodule, and the quotient to the next comodule
	matrix_mem<BP> inj;
	matrix_file<BP> indj, qut, new_map;
	//the matrices for the resolutions
	std::vector<matrix_mem<Z3>> mapses;
	matrix_file<Z3> mm;
	
	//the curtis table for the resolutions
	std::vector<curtis_table_mem<F3>> ResolutionTables;
	curtis_table_mem<F3> ctable;
	
	//the generators for a resolution
	std::vector<std::vector<int>> gens;
	
	//a comdule which is initialization to the trivial comodule of rank 1
	BPComodInit comod;
	
	//the operations for multiplicative structures
	multiplication multp;
	
	//load the curtis table data
	void loadResolutionTables(string table_data);
	
	//load the data for generators
	void load_gens(string gens_data);
	
	//Write a self-describing record of this run beside its tables, so that
	//downstream tools do not have to guess which ring it was done over.
	//Called from the constructor, so every run has one.
	void save_run_info();
	
	//the constructor. hgt defaults to 0, so every existing caller keeps the
	//classical BP_* computation unchanged.
	BPInit(int max_deg, int resolution_length, string etaL_data, string delta_data, string R2L_data, string dirname, int hgt = 0);
	
	//do resolutions
	void resolve();
	
	//construct resolution
	void resolution();
	
	//make algebraic Novikov table
	void make_algNov();
	//make Bockstein table
	void make_Boc();
	
	//make multiplication table
	void mult_table(BPBP const &, int, string);
	void mult_table();
	//make multiplication table for top theta on the Moore spectrum
	void mult_theta(int);
};
