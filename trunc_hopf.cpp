//trunc_hopf.cpp
#include "trunc_hopf.h"

//constructor. maxdeg is the total internal-degree range the resolution will be computed
//over (same role as Steenrod_Op's maxdeg / mr_st's CLI max_degree argument) -- it must be
//set at least as large as the highest degree the requested resolution_length will reach.
TruncHopf_Op::TruncHopf_Op(int maxdeg, int prime) : Fp_Op(prime), Gamma_opers(this), FpMod_opers(this), GammaMod_opers(&Gamma_opers){
	this->maxDeg = maxdeg;

	Hopf_Algebroid<Fp,Gamma>::ringOper = this;
	algebroidRingOper = &Gamma_opers;

	moduleOper = &FpMod_opers;
	algebroidModuleOper = &GammaMod_opers;
}

//Gamma's basis element x^n is represented as the exponent n itself (single variable, no
//packing needed), so the basis index (matrix_index) and the exponent coincide.
vectors<matrix_index, Fp> TruncHopf_Op::algebroid2vector(const Gamma& x, int shift){
	std::function<matrix_index(exponent)> rd = [shift](exponent e){
		return (matrix_index)e + shift; };
	return Gamma_opers.re_index(rd, x);
}

Gamma TruncHopf_Op::vector2algebroid(const vectors<matrix_index, Fp>& v){
	std::function<exponent(matrix_index)> rd = [](matrix_index n){
		return (exponent) n; };
	return FpMod_opers.re_index(rd, v);
}

Gamma TruncHopf_Op::etaL(const Fp& x){
	return Gamma_opers.monomial(0,x); }

Gamma TruncHopf_Op::etaR(const Fp& x){
	return Gamma_opers.monomial(0,x); }

//binomial coefficient mod p, 0<=i<=n<GammaDim (Pascal's triangle, computed once per call --
//GammaDim is tiny so this is not a performance concern)
Fp TruncHopf_Op::binom(int n, int i){
	if(i<0 || i>n) return this->zero();
	std::vector<std::vector<int>> table(n+1, std::vector<int>(n+1,0));
	for(int r=0; r<=n; ++r){
		table[r][0] = 1;
		for(int c=1; c<=r; ++c)
			table[r][c] = table[r-1][c-1] + (c<=r-1 ? table[r-1][c] : 0);
	}
	return this->unit(table[n][i]);
}

//reduced comultiplication of x^n: Delta(x^n) = sum_{i=0}^n binom(n,i) x^i (x) x^{n-i},
//extended multiplicatively from Delta(x) = 1(x)x + x(x)1
vectors<matrix_index, Gamma> TruncHopf_Op::delta(matrix_index n){
	vectors<matrix_index, Gamma> result;
	for(int i=0; i<=(int)n; ++i){
		Fp c = binom((int)n,i);
		if(!this->isZero(c))
			result.push({(matrix_index)i, Gamma_opers.monomial(n-i,c)});
	}
	return result;
}

//the number of basis elements (1,x,...,x^{GammaDim-1}) of degree strictly below n
unsigned TruncHopf_Op::ranksBelowDeg(unsigned n){
	return n<GammaDim ? n : GammaDim; }

//initialization: Gamma's own degree grading is just the exponent itself
void TruncHopf_Op::initialize(){
	std::function<int(matrix_index)> cofree_degs = [](matrix_index n){
		return (int) n; };
	this->init_cofree_data(cofree_degs);
}

//the underlying degree
template<>
int FreeTruncCoMod::underlyingDeg(int i){
	return i; }

//add degrees
template<>
int FreeTruncCoMod::add_degree(int a, int b){
	return a+b; }

//output the generators of a resolution: number of generators per (internal degree, step)
string TruncHopf_Op::output_resolution(std::iostream &gens_file, int length){
	std::vector<FreeTruncCoMod> gens(length);
	for(int i=0; i<length-1; ++i){
		int32_t M_rank;
		gens_file.read((char*)&M_rank, 4);
		gens[i].load(gens_file);
	}

	string result = "\n";

	for(int i=length-1; i>=0; --i)
		result += gens[i].generators.output() + "\n";

	result += "\n";

	std::vector<std::vector<int>> dims(maxDeg, std::vector<int>(length));

	for(int i=0; i<length; ++i)
		for(auto d: gens[i].generators.degree)
			dims[d-i][i]++;

	for(int i=length-1; i>=0; --i){
		for(int j=0; j<(int)maxDeg; ++j)
			result += std::to_string(dims[j][i]) + "\t";
		result += "\n";
	}

	return result;
}

string TruncHopf_Op::output_resolution(string gens_file, int length){
	std::fstream gf(gens_file, std::ios::in | std::ios::binary);
	return output_resolution(gf,length);
}
