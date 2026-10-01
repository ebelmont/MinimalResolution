//BP_init.cpp
#include"BP_init.h"

//initialize a generic comodule
BPComodInit::BPComodInit(string filename) : BPCoMod_generic(&coaction_matrix), coaction_matrix(filename){
}

//the constructor
//
//Note the base ring's arithmetic is chosen here, before anything else is
//built: at height 0 it is Z_(p) as always, above it F_p, because
//BP_*/I_n = F_p[v_n, v_{n+1}, ...] has characteristic p. Everything that
//consumes the base ring does so through a RingOp<Z3>*, so this one choice
//propagates through BP, BPBP, and every module and matrix over them.
BPInit::BPInit(int max_deg, int res_length, string etaL_data, string delta_data, string R2L_data, string dirname, int hgt) : Z3_oper_store(hgt > 0 ? (Z3_Op*)new Z3_mod_p_Op() : new Z3_Op()), Z3_oper(*Z3_oper_store), BP_oper(max_deg, &Z3_oper, &etaL_matrix, &delta_matrix, &R2L_matrix), F3Mod_opers(&Z3_oper.F3_opers), etaL_matrix(dirname + "etaL_matrix"), R2L_matrix(dirname + "R2L_matrix"), delta_matrix(dirname + "delta_matrix"), indj(dirname + "indj"), qut(dirname + "qut"), new_map(dirname + "new_map"), mm(dirname + "mm_matrix"), comod(dirname + "comodule_matrix"), multp(&BP_oper){
	max_degree = max_deg;
	resolution_length = res_length;
	director = dirname;
	height = hgt;
	
	//initialize matric operators
	matrix<BP>::moduleOper = &BP_oper.BPMod_opers;
	matrix<BPBP>::moduleOper = &BP_oper.BPBPMod_opers;
	matrix<Z3>::moduleOper = &BP_oper.Z3Mod_opers;

	//The height has to be set BEFORE the structure tables are read: they are
	//reduced mod I_n as they stream in (matrix_file's update_all is
	//unimplemented, so there is no reducing them after the fact).
	BP_oper.set_height(height);
	
	//initialize the structure data
	BP_oper.initialize(etaL_data, delta_data, R2L_data);
	
	//initialize the comod to a trivial one with one generator at degree 0.
	//At height n this is BP_*/I_n as a comodule over BP_*BP/I_n -- i.e.
	//already the right thing for S/p, S/(p,v_1), ... -- exactly as it is
	//BP_* itself (the sphere) at height 0.
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
	
	//record what ring this run was done over, for whoever reads the tables
	save_run_info();
}

//Write a machine-readable record of this run beside its tables.
//
//WHY THIS EXISTS. The output tables do not say which ring they were computed
//over, and the answer changes what they MEAN: Ext^0 is torsion free over
//BP_*, but an F_p-vector space over BP_*/I_n, so the same list of class names
//denotes Z_(p)'s in one case and Z/p's in the other. Without this file a
//consumer has to guess -- from the output filenames, or from whether any
//class name happens to mention v_0 -- and a wrong guess is silently wrong,
//not an error. anss_chart.py reads this; see docs/CHARTS.md.
//
//Format: one "key value" per line, '#' starts a comment.
void BPInit::save_run_info(){
	std::fstream info(director + "run_info.txt", std::ios::out);
	if(!info.is_open()){
		std::cerr << "warning: could not write " << director << "run_info.txt\n" << std::flush;
		return;
	}
	info << "# what this run computed. Written by BPInit; read by anss_chart.py.\n";
	info << "height " << height << "\n";
	//The characteristic of the base ring -- 0 for BP_*, p for BP_*/I_n with
	//n >= 1. This is the fact consumers actually need, so publish it rather
	//than making each of them re-derive it from the height.
	info << "characteristic " << (height > 0 ? Z3_oper.F3_opers.prime : 0) << "\n";
	info << "max_degree " << max_degree << "\n";
	info << "resolution_length " << resolution_length << "\n";
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
	
	//Multiplication by the bottom generator of the base ring's maximal
	//invariant ideal: p = v_0 over BP_*, v_n over BP_*/I_n (where p is 0, so
	//the v_0 table would be uniformly zero and tell you nothing).
	if(height > 0){
		auto et2 = multp.vn_extension(resolution_length,Complex,AANtables,height,resolution_length);
		std::fstream anf(director + "AANSS_a" + std::to_string(height) + ".txt", std::ios::out);
		anf << multp.output_multiplication_table(et2,0,resolution_length-1);
	}
	else{
		auto et2 = multp.three_extension(resolution_length,Complex,AANtables,resolution_length);
		std::fstream a0f(director + "AANSS_a0.txt", std::ios::out);
		a0f << multp.output_multiplication_table(et2,0,resolution_length-1);
	}
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
	
	//as in make_algNov: p over BP_*, v_n over BP_*/I_n
	if(height > 0){
		auto et2 = multp.vn_extension(resolution_length,Complex,Btables,height,resolution_length);
		std::fstream anf(director + "BocSS_a" + std::to_string(height) + ".txt", std::ios::out);
		anf << multp.output_multiplication_table1(et2,0,resolution_length-1);
	}
	else{
		auto et2 = multp.three_extension(resolution_length,Complex,Btables,resolution_length);
		std::fstream a0f(director + "BocSS_a0.txt", std::ios::out);
		a0f << multp.output_multiplication_table1(et2,0,resolution_length-1);
	}
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
