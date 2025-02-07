//BP_init.cpp
#include"BP_init.h"

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
		
	BP_oper.pre_resolution_modeled(comod, director + "maps", director + "gens", resolution_length, director+"tables", &ctable, gens, tfm, &inj, &qut, &indj, &new_map, director + "back");
	
	//combining generator files
	BP_oper.gens_file_combiner(director + "gens", resolution_length, comod);
}

//construct resolution
void BPInit::resolution(){
	//construct the resolution
	std::fstream outfile("maps.txt", std::ios::out);
	BP_oper.resolution(director + "maps", director + "gens", director + "res" , resolution_length, comod, &inj, &qut, &indj, NULL);
	
	//set the matrices
	mapses.resize(resolution_length + 2);
	std::function <matrix<Z3>*(int)> mst = [this](int i){
		return &mapses[i]; };
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

//make algebraic Novikov table
void BPInit::make_algNov(){
	Complex.load_matrix(resolution_length,director + "gens", director+"cpx");
	
	AANtables.table_of_complex(Complex,resolution_length);
	AANtables.save(director + "AANSS_table_binary");
	std::fstream atb(director + "AANSS_table.txt", std::ios::out);
	atb << AANtables.output_tables();
	
	auto et2 = multp.three_extension(resolution_length,Complex,AANtables,resolution_length);
	std::fstream a0f(director + "AANSS_a0.txt", std::ios::out);
	a0f << multp.output_multiplication_table(et2,0,resolution_length-1);
}

//make Bockstein table
void BPInit::make_Boc(){
	Complex.load_matrix(resolution_length,director + "gens", director+"cpx");
	
	Btables.table_of_complex(Complex,resolution_length);
	Btables.save(director + "BocSS_table_binary");
	std::fstream btb(director + "BocSS_table.txt", std::ios::out);
	btb << Btables.output_tables();
	
	auto ba = Btables.Bname2Anames(resolution_length, AANtables, resolution_length);
	std::fstream b2a(director + "B2A_table.txt", std::ios::out);
	b2a << multp.output_multiplication_table(ba,0,resolution_length-1);
	
	auto et2 = multp.three_extension(resolution_length,Complex,Btables,resolution_length);
	std::fstream a0f(director + "BocSS_a0.txt", std::ios::out);
	a0f << multp.output_multiplication_table1(et2,0,resolution_length-1);
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
	//compute the table for the multiplication on algebraic Novikov spectral sequence
	auto multable = multp.mult_extension(&inj, max_degree-deg, resolution_length, genst, director + "res", Complex, AANtables, &indj, &mm, resolution_length);
	//output the table
	std::fstream file(director + "AANSS_" + filename, std::ios::out);
	file << multp.output_multiplication_table(multable, 1, resolution_length-2);
	//compute the table for the multiplication on Bockstein spectral sequence
	auto Bmultable = multp.mult_extension1(&inj, max_degree-deg, resolution_length, genst, director + "res", Complex, Btables, &indj, &mm, 1, true);
	//output the table
	std::fstream fileB(director + "BocSS_" + filename, std::ios::out);
	fileB << multp.output_multiplication_table(Bmultable, 1, resolution_length-2);
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
		//compute the table for the multiplication on Bockstein spectral sequence
		auto Bmultable = multp.mult_extension(&inj, max_degree-deg[i-2], resolution_length, genst, director + "res", Complex, Btables, &indj, &mm, 1, true);
		//output the table
		std::fstream fileB(director + "BocSS_theta" + std::to_string(i) + ".txt", std::ios::out);
		fileB << multp.output_multiplication_table(Bmultable, 1, resolution_length+1); // EB: the last argument here is 1 + the last homological degree where theta products are recorded in the output files ...theta2.txt, etc.
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
