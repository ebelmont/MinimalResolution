//mr_trunc.cpp
//computes the minimal cofree resolution of F_p as the trivial comodule over the Hopf
//algebroid (F_p, Gamma), Gamma = F_p[x]/(x^GammaDim - 1), |x|=1, Delta(x)=1(x)x+x(x)1
//extended multiplicatively. See docs/GENERAL_COMODULES.md-adjacent notes for context;
//unlike the BP_*BP case, F_p is already a field so pre_resolution_tab applies directly,
//with no mod-I reduction / lifting needed.
#include"trunc_init.h"

int main(int argc, char** argv){
	string filename = argv[1];
	filename += "_";

	int maxdeg = std::atoi(argv[1]);
	int length = std::atoi(argv[2]);
	int prime = argc>3 ? std::atoi(argv[3]) : 7;

	TruncInit tr(prime, maxdeg, length);

	tr.resolve(filename);
	tr.saveResolutionTables(filename + "ResTables");
	tr.save_gens(filename + "gens_data");

	std::cout << tr.trunc_oper.output_resolution(filename+"gens",length);

	return 0;
}
