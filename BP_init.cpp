//BP_init.cpp
#include"BP_init.h"
#include<limits>
#include<algorithm>

//initialize a generic comodule
BPComodInit::BPComodInit(string filename) : BPCoMod_generic(&coaction_matrix), coaction_matrix(filename){
}

//the constructor
BPInit::BPInit(int max_deg, int res_length, string etaL_data, string delta_data, string R2L_data, string dirname) : BP_oper(max_deg, &Z3_oper, &etaL_matrix, &delta_matrix, &R2L_matrix), F3Mod_opers(&Z3_oper.F3_opers), etaL_matrix(dirname + "etaL_matrix"), R2L_matrix(dirname + "R2L_matrix"), delta_matrix(dirname + "delta_matrix"), indj(dirname + "indj"), qut(dirname + "qut"), new_map(dirname + "new_map"), mm(dirname + "mm_matrix"), comod(dirname + "comodule_matrix"), multp(&BP_oper){
	max_degree = max_deg;
	resolution_length = res_length;
	director = dirname;
	
	//initialize matric operators
	matrix<BP>::moduleOper = &BP_oper.BPMod_opers;
	matrix<BPBP>::moduleOper = &BP_oper.BPBPMod_opers;
	matrix<Z3>::moduleOper = &BP_oper.Z3Mod_opers;

	//initialize the structure data
	BP_oper.initialize(etaL_data, delta_data, R2L_data);
	
	//initialize the comod to a trivial one with one generator at degree 0
	BP_oper.set_to_trivial(comod,0);
	
	//initialize the curtis tables
	curtis_table<F3>::ModOper = &F3Mod_opers;
	ResolutionTables.resize(resolution_length+2);
	
	//set the complex of primitives
	primitive_data::set_oper(&BP_oper);
	
	//set the algebraic Novikov table
	algNov_table::set_op(&BP_oper);
	AAN_table.resize(resolution_length+1);
	AANtables.set_table(&AAN_table);
	
	//set the Bockstein table
	B_table.resize(res_length+1);
	Btables.set_table(&B_table);
}

//do resolutions
void BPInit::resolve(){
	std::function<curtis_table<F3>*(int)> tables = [this](int i){
		return &ResolutionTables[i]; };
	
	std::function<vectors<matrix_index,BP>(const vectors<matrix_index,Fp>&)> tfm = [this] (const vectors<matrix_index,Fp>& v){
		return BP_oper.lift(v); };

	//resume from the highest-numbered back-up file already on disk (each
	//back<N> holds the resolved comodule's state *after* degree N was
	//completed), instead of always starting the resolution over at degree 0
	std::string back_prefix = director + "back";
	int res_start = -1;
	for(int i=0; i<=resolution_length; ++i){
		std::ifstream back_file(back_prefix + std::to_string(i), std::ios::in | std::ios::binary);
		if(!back_file.is_open()) break;
		res_start = i;
	}

	BP_oper.pre_resolution_modeled(comod, director + "maps", director + "gens", resolution_length, director+"tables", &ctable, gens, tfm, &inj, &qut, &indj, &new_map, back_prefix, res_start);

	//combining generator files
	BP_oper.gens_file_combiner(director + "gens", resolution_length, comod);
}

//construct resolution
void BPInit::resolution(){
	//Complex.Maps must hold one allocated matrix<Z3> per resolution degree
	//before make_algNov()'s Complex.load_matrix() call, whether or not the
	//expensive rebuild below runs: load_matrix() only fills in already-
	//existing Maps[i] objects, it doesn't construct them. Cheap (no disk
	//I/O), so just always do it up front.
	mapses.resize(resolution_length + 2);
	std::function <matrix<Z3>*(int)> mst = [this](int i){
		return &mapses[i]; };
	Complex.Maps.resize(resolution_length + 1);
	for(int i=0; i<=resolution_length; ++i)
		Complex.Maps[i] = mst(i);

	//if a complex-of-primitives matrix is already on disk from a prior
	//completed run, don't rebuild it: this step's only other lasting
	//output is `cpx` (the `res`/`maps`/`gens` files it reads are not read
	//by anything downstream -- make_algNov()/make_Boc() each reload the
	//complex fresh from `cpx` themselves via their own
	//Complex.load_matrix() calls), so an existing `cpx` means this whole,
	//expensive step can be skipped entirely
	std::ifstream cpx_check(director + "cpx", std::ios::in | std::ios::binary);
	if(cpx_check.is_open() && cpx_check.peek() != std::ifstream::traits_type::eof())
		return;
	cpx_check.close();

	//construct the resolution
	std::fstream outfile("maps.txt", std::ios::out);
	BP_oper.resolution(director + "maps", director + "gens", director + "res" , resolution_length, comod, &inj, &qut, &indj, NULL);

	//construct the complex of primitives
	Complex.load(resolution_length,director + "gens",director + "res",&inj,&indj,mst);
	Complex.save_matrix(director + "cpx");
}

//load the curtis table data
void BPInit::loadResolutionTables(string table_data){
	std::fstream tables(table_data, std::ios::in | std::ios::binary);
	std::cout << tables.is_open() << std::flush;
	for(unsigned i=0; i<ResolutionTables.size(); ++i){
		ResolutionTables[i].load(tables);
		std::cout << ResolutionTables[i].output();
	}
}

//load the data for generators
void BPInit::load_gens(string gens_data){
	std::fstream genfile(gens_data, std::ios::in | std::ios::binary);
	std::cout << genfile.is_open() << std::flush;
	int32_t sz;
	genfile.read((char*)&sz, 4);
	gens.resize(sz);
	for(int i=0; i<sz; ++i){
		int32_t ss;
		genfile.read((char*)&ss, 4);
		gens[i].resize(ss);
		for(int j=0; j<ss; ++j)
			genfile.read((char*)&gens[i][j], 4);
	}
}

//No cap on the algebraic Novikov filtration f (total power of v0,v1,v2,...
//in a class's leading monomial): every class with s < resolution_length and
//t <= max_degree is searched and printed. The cap used to be
//resolution_length - s, i.e. on s+f rather than on s, which silently dropped
//classes with t and s in range but large f (e.g. v1^40[1-0] at
//(stem,s)=(163,1) with resolution_length 40). t already bounds the
//v1,v2,... part and Z3's precision (Z3_Op::MAX_VALUATION) bounds the v0 part.
//The filtration-limit parameters themselves stay because the Bockstein
//product tables use them for real (limit 0, in mult_table/mult_theta).
static const int NO_FILTRATION_LIMIT = std::numeric_limits<int>::max();

//make algebraic Novikov table
void BPInit::make_algNov(){
	Complex.load_matrix(resolution_length,director + "gens", director+"cpx");

	//pric bounds total algebraic-Novikov filtration search depth (v0- AND
	//v1,v2,...-power combined) -- unrelated to Z3's own v0-only precision
	//ceiling, so it must NOT be clamped to Z3_Op::MAX_VALUATION (that
	//wrongly truncated genuine high-v1-power towers with zero v0 content;
	//see NOTES.md "pric conflated total filtration with v0 precision").
	//The actual v0-precision guard lives in cancel_scalar (SS.h), which
	//checks the real v0 valuation directly and hard-aborts if it would
	//exceed Z3_Op::MAX_VALUATION.
	int pric = NO_FILTRATION_LIMIT;

	AANtables.table_of_complex(Complex,pric,director + "AANback");
	AANtables.save(director + "AANSS_table_binary");
	std::fstream atb(director + "AANSS_table.txt", std::ios::out);
	atb << AANtables.output_tables(pric);

	auto et2 = multp.three_extension(resolution_length,Complex,AANtables,pric);
	std::fstream a0f(director + "AANSS_a0.txt", std::ios::out);
	a0f << multp.output_multiplication_table(et2,0,pric);
}

//make Bockstein table
void BPInit::make_Boc(){
	Complex.load_matrix(resolution_length,director + "gens", director+"cpx");

	//see make_algNov's comment on pric vs resolution_length
	int pric = NO_FILTRATION_LIMIT;

	Btables.table_of_complex(Complex,pric);
	Btables.save(director + "BocSS_table_binary");
	std::fstream btb(director + "BocSS_table.txt", std::ios::out);
	btb << Btables.output_tables(pric);

	auto ba = Btables.Bname2Anames(resolution_length, AANtables, pric);
	std::fstream b2a(director + "B2A_table.txt", std::ios::out);
	b2a << multp.output_multiplication_table(ba,0,pric);

	auto et2 = multp.three_extension(resolution_length,Complex,Btables,pric);
	std::fstream a0f(director + "BocSS_a0.txt", std::ios::out);
	a0f << multp.output_multiplication_table1(et2,0,pric);
}

//make the algebraic Novikov multiplication table by a given element in BPBP
void BPInit::mult_table(BPBP const &mul, int deg, string filename){
	//load the complex
	auto genst = BPComplex::get_generator(resolution_length, director + "gens");
	//load the complex of primitives
	Complex.load_matrix(resolution_length, director +"gens", director + "cpx");
	//load the algebraic Novikov table
	AANtables.load(director + "AANSS_table_binary");
	Btables.load(director + "BocSS_table_binary");
	//make the multiplication table on BPBP
	multp.make_eta_R_multiplier(mul, &inj, deg);
	//compute the table for the multiplication on algebraic Novikov spectral sequence.
	//see make_algNov's comment on pric vs resolution_length -- not the
	//Z3 precision ceiling, which cancel_scalar guards independently.
	int pric = NO_FILTRATION_LIMIT;
	auto multable = multp.mult_extension(&inj, max_degree-deg, resolution_length, genst, director + "res", Complex, AANtables, &indj, &mm, pric);
	//output the table
	std::fstream file(director + "AANSS_" + filename, std::ios::out);
	file << multp.output_multiplication_table(multable, 1, pric);
	//compute the table for the multiplication on Bockstein spectral sequence.
	//Limit 0 (inclusive, like every pric): name products through v0-power 0
	//only, i.e. the action on the Bockstein E_1 = Ext(BP_*/3).
	auto Bmultable = multp.mult_extension1(&inj, max_degree-deg, resolution_length, genst, director + "res", Complex, Btables, &indj, &mm, 0, true);
	//output the table
	std::fstream fileB(director + "BocSS_" + filename, std::ios::out);
	fileB << multp.output_multiplication_table(Bmultable, 1, NO_FILTRATION_LIMIT);
}

//make multiplication table for top theta on the Moore spectrum
void BPInit::mult_theta(int resolution_length){
	//load the complex
	auto genst = BPComplex::get_generator(resolution_length, director + "gens");
	//load the complex of primitives
	Complex.load_matrix(resolution_length, director +"gens", director + "cpx");

	//load the algebraic Novikov table
	AANtables.load(director + "AANSS_table_binary");
	Btables.load(director + "BocSS_table_binary");
	
	//compute the top thetas
	auto theta = BP_oper.thetas();
	
	//degree of theta
	int deg[maxVar+1] = {12,28,36,60,76,84};
	for(int i=2; i<=7; ++i){
		//make the multiplication table on BPBP
		multp.make_eta_R_multiplier(theta[i], &inj, deg[i-2]);
		//compute the table for the multiplication on Bockstein spectral sequence.
		//Limit 0 (inclusive): keep only the v0^0 terms of theta*x. beta*x =
		//delta(theta*x) depends only on those -- a v0^k term (k>=1) differs
		//from zero by 3*(chain), whose delta is a boundary -- so higher terms
		//would only add work downstream and then cancel.
		auto Bmultable = multp.mult_extension(&inj, max_degree-deg[i-2], resolution_length, genst, director + "res", Complex, Btables, &indj, &mm, 0, true);
		//output the table
		std::fstream fileB(director + "BocSS_theta" + std::to_string(i) + ".txt", std::ios::out);
		fileB << multp.output_multiplication_table(Bmultable, 1, NO_FILTRATION_LIMIT); // EB: the last argument here is 1 + the last homological degree where theta products are recorded in the output files ...theta2.txt, etc. [It is actually a filtration cutoff on the multiplied class; see NO_FILTRATION_LIMIT.]
	}
}

//make the algebraic Novikov multiplication table by a given element in BPBP
void BPInit::mult_table(){
	//make the h0 multiplication table
	BPBP h0 = BP_oper.h0();
	// h0 is defined as (eta_R(v1) - eta_L(v1))/p
	BPBP h02 = BP_oper.BPBP_opers.multiply(h0, h0);
	mult_table(h0, 4, "h0.txt");
//	mult_table(h02, 8, "h0squared.txt");
}
