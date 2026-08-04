//mr_BP_comod.cpp
//
//Computes the algebraic Novikov E2 page of a BP_*BP-comodule chosen by name
//on the command line. The comodules themselves live in comodules.cpp -- add
//new ones there, not here.
//
//    ./mr_BP_comod <halfT> <resolution_length> [comodule]
//    ./mr_BP_comod --list
//
//<comodule> defaults to "sphere" (the trivial comodule BP_*), so the default
//invocation computes the same thing the shipped mr_BP does.
//
//Needs a matching `BPtab <halfT>` run first, exactly as mr_BP does -- that
//supplies the (BP_*, BP_*BP) structure maps, which describe the Hopf
//algebroid itself and so are the same whichever comodule is being resolved.
//See docs/GENERAL_COMODULES.md and docs/BUILD_AND_RUN.md.
//
//HOW IT WORKS: a resolution cannot be searched for directly over BP_* (see
//BP_generic_init.h for why), so this reduces the chosen comodule mod
//I = (p,v_1,v_2,...) to a comodule over the FIELD P = BP_*BP/I, resolves
//that with the classical Steenrod-style machinery, then lifts that
//resolution back to BP_*BP. You describe your comodule once, in terms of
//BP_*BP; reduce_coaction_rows_mod_I derives the P-side reduction.
//
//Output files are prefixed <halfT>_<comodule>BP... (final resolution,
//algNov/Boc tables) and <halfT>_<comodule>P... (the mod-I model resolution),
//so runs for different comodules -- and a plain mr_st/mr_BP run -- coexist
//in one directory without clobbering each other.

#include"BP_generic_init.h"
#include"Steenrod_generic_init.h"
#include"BP_mod_I.h"
#include"comodules.h"

static void usage(const char *prog){
	std::cout << "usage: " << prog << " <halfT> <resolution_length> [comodule]\n"
	          << "       " << prog << " --list\n\n"
	          << "available comodules:\n" << list_comodules()
	          << "\nRun `BPtab <halfT>` first, as for mr_BP.\n" << std::flush;
}

int main(int argc, char** argv){
	if(argc >= 2 && (string(argv[1]) == "--list" || string(argv[1]) == "-l")){
		std::cout << "available comodules:\n" << list_comodules() << std::flush;
		return 0;
	}
	if(argc >= 2 && (string(argv[1]) == "--help" || string(argv[1]) == "-h")){
		usage(argv[0]);
		return 0;
	}
	if(argc < 3){
		usage(argv[0]);
		return 1;
	}

	int max_degree = std::atoi(argv[1]);
	int resolution_length = std::atoi(argv[2]);
	string comod_name = argc > 3 ? string(argv[3]) : default_comodule_name();

	if(max_degree <= 0 || resolution_length <= 0){
		std::cerr << "error: <halfT> and <resolution_length> must be positive integers\n\n" << std::flush;
		usage(argv[0]);
		return 1;
	}

	const ComoduleSpec *spec = find_comodule(comod_name);
	if(spec == NULL){
		std::cerr << "error: unknown comodule \"" << comod_name << "\"\n\n"
		          << "available comodules:\n" << list_comodules() << std::flush;
		return 1;
	}

	std::cout << "comodule: " << spec->name << " -- " << spec->description << "\n" << std::flush;

	string filename0  = string(argv[1]) + "_";
	string bp_dir     = filename0 + spec->name + "BP";   //final (BP-side) output prefix
	string model_dir  = filename0 + spec->name + "P";    //mod-I model output prefix

	//construct the BP-side driver first: this loads the structure maps a
	//BPtab run already produced, so the builder below can use BP_oper.h0()
	//and friends to build BP_*BP elements
	BPGenericInit BPoper(max_degree, resolution_length, filename0+"etaL", filename0+"R2L", filename0+"delta", bp_dir);

	int rank;
	std::vector<int> degree;
	std::function<vectors<matrix_index,BPBP>(int)> coaction_rows;
	spec->build(BPoper.BP_oper, rank, degree, coaction_rows);

	std::cout << "rank " << rank << ", generator degrees:";
	for(auto d : degree) std::cout << " " << d;
	std::cout << "\n" << std::flush;

	BPoper.set_comodule(rank, degree, coaction_rows);

	//phase 1: resolve the comodule's reduction mod I over the field-based
	//Hopf algebra P = BP_*BP/I. The model needs one more step than the
	//BP-side length, mirroring the mr_st/mr_BP convention.
	SteenrodGenericInit stOper(3, max_degree, resolution_length+1, model_dir+"steenrod_coaction.data");
	stOper.set_comodule(rank, degree, reduce_coaction_rows_mod_I(coaction_rows));
	stOper.resolve(model_dir);
	stOper.save_gens(model_dir + "gens_data");

	//phase 2: lift that model resolution to a BP_*BP-comodule resolution
	std::cout << "starting resolution..." << std::flush;
	BPoper.resolve(model_dir + "gens_data", model_dir + "BPtables");

	//post-processing, unchanged from mr_BP -- these only read back the
	//resolution's own output files
	BPoper.resolution();
	BPoper.make_algNov();
	BPoper.make_Boc();

	std::cout << "\ndone. output written with prefix " << bp_dir << "\n" << std::flush;
	return 0;
}
