//trunc_hopf.h
#pragma once
#include "Fp.h"
#include "polynomial.h"
#include "matrices.h"
#include "hopf_algebroid.h"

//Gamma = F_p[x]/(x^GammaDim - 1), |x| = 1, i.e. F_p[x] truncated above degree GammaDim-1.
//Represented as an ordinary single-variable polynomial<Fp>; the truncation is enforced by
//ranksBelowDeg (which never returns more than GammaDim) and by delta (which only ever asks
//for basis elements 0..GammaDim-1). Ordinary (untruncated) multiplication in Gamma_opers is
//only ever invoked by the resolution machinery with a degree-0 (scalar) factor -- see
//Hopf_Algebroid::right_scalor_mult -- so no special truncated multiply is needed.
constexpr unsigned GammaDim = 5;

//the algebroid Gamma
typedef polynomial<Fp> Gamma;

class TruncHopf_Op;

//the class of the Hopf algebroid (F_p, Gamma)
class TruncHopf_Op : virtual public Hopf_Algebroid<Fp,Gamma>, virtual public Fp_Op{
public:
	//ring operations on Gamma
	PolynomialOp_Para<Fp> Gamma_opers;

	//operations on Fp-modules
	ModuleOp<matrix_index,Fp> FpMod_opers;

	//operations on Gamma-modules
	ModuleOp<matrix_index,Gamma> GammaMod_opers;

	//the constructor
	TruncHopf_Op(int maxdeg, int prime);

	//transform an element of Gamma to a vector over Fp
	vectors<matrix_index, Fp> algebroid2vector(const Gamma&,int shift);
	//the inverse of the transformation
	Gamma vector2algebroid(const vectors<matrix_index, Fp>&);

	//the left unit
	Gamma etaL(const Fp&);
	//the right unit
	Gamma etaR(const Fp&);

	//the co-multiplication on the n-th basis element x^n, i.e. reduced Delta(x^n) = sum_i binom(n,i) x^i (x) x^(n-i)
	vectors<matrix_index, Gamma> delta(matrix_index n);

	//the number of basis elements (1,x,...,x^{GammaDim-1}) of degree below a given bound
	unsigned ranksBelowDeg(unsigned);

	//the initialization
	void initialize();

	//output the number of generators in each degree of a resolution
	string output_resolution(std::iostream &gens_file, int length);
	string output_resolution(string gens_file, int length);

private:
	//binomial coefficient mod p, for 0<=i<=n<GammaDim
	Fp binom(int n, int i);
};

//the class of cofree comodules over Gamma
typedef cofree_comodule<Gamma,int> FreeTruncCoMod;

//a comodule over Gamma
typedef comodule_generic<Gamma,int> TruncCoMod_generic;
