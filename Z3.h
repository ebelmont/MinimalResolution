//Z3.h
#pragma once
#include"algebra.h"
#include"Fp.h"
#include <cstdint>

//we use 64-bit unsigned integer to denote 3-adic numbers
typedef uint64_t Z3;

//the residue field
typedef Fp F3;

class Z3_Op : public virtual RingOp<Z3>{
public:
	//the constructor
	Z3_Op();
	
	//the charateristic of the residue field
	int prime();
	
	//the operations on F3
	Fp_Op F3_opers;
	
	//the valuation
	unsigned valuation(Z3);
	
	//n-th power of p
	Z3 power_p(int n);
	
	// divide by p^n
	Z3 divide(const Z3&, int n); 
	
	//lift from the residue field
	Z3 lift(F3);
    
	//addition
	Z3 add(const Z3&, const Z3&);
	Z3 add(Z3&&, Z3&&);
    
	//multiplication
	Z3 multiply(const Z3&, const Z3&);
    
	//the unit
	Z3 unit(int);
    
	//check if it is zero, or rather divisible by 2^64
	bool isZero(const Z3&);
	
	//the zero element
	Z3 zero();
	
	//negation
	Z3 minus(const Z3&);
    
	//IO operations
	string output(Z3);
	void save(const Z3&, std::iostream&);
	Z3 load(std::iostream&);
    
	//check if invertible
	virtual bool invertible(const Z3&);
	//inverse of an invertible element
	virtual Z3 inverse(const Z3&);
	
	//the canonical representative of an element. The identity here; the
	//mod-p subclass below reduces. Anything that has to decide whether a
	//coefficient read off a file or a table is really zero should push it
	//through this first -- see BP_Op's reduce_*_mod_I helpers.
	virtual Z3 normalize(const Z3&);
	
	//this class is used polymorphically (BPInit owns one through a
	//Z3_Op pointer), so it needs a virtual destructor
	virtual ~Z3_Op(){}
};

//Z3_Op with the arithmetic of F_p = Z_(p)/(p) instead of that of Z_(p).
//
//The BP_*/I_n pipeline needs a base ring of characteristic p, because
//BP_*/I_n = F_p[v_n, v_{n+1}, ...] for n >= 1. Rather than introduce a second
//coefficient TYPE, this keeps the Z3 representation and changes only the
//arithmetic, so that every element is normalised to {0,...,p-1}.
//
//Why not just use polynomial<Fp>: that type is already taken -- it is P, the
//dual-Steenrod-side algebra -- and matrix<P>::moduleOper is a static that
//both halves of mr_BP_comod share inside one process, so a second use of the
//type would clobber it.
//
//Substituting this for Z3_Op when BP_Op is constructed is enough to turn
//BP_* into BP_*/p everywhere: BP and BPBP arithmetic both run through
//ringOper, and so does every module and matrix operation above them. What it
//does NOT do is kill v_1,...,v_{n-1}; that is a statement about exponents,
//and lives in BP_Op::reduce_*_mod_I.
//
//The inherited non-virtual helpers still behave correctly on normalised
//values: valuation() is 0 on a nonzero residue and max_val on 0 (it is built
//from the virtual isZero/invertible), and power_p(n) is 0 for n > 0 (it is
//built from unit(p)). divide() is NOT meaningful here -- see BP_Op::h0(),
//which takes a different route in characteristic p.
class Z3_mod_p_Op : public Z3_Op{
public:
	//addition
	Z3 add(const Z3&, const Z3&);
	Z3 add(Z3&&, Z3&&);
	
	//multiplication
	Z3 multiply(const Z3&, const Z3&);
	
	//the unit map from the integers
	Z3 unit(int);
	
	//check if it is zero
	bool isZero(const Z3&);
	
	//negation
	Z3 minus(const Z3&);
	
	//invertibility: every nonzero residue is a unit
	bool invertible(const Z3&);
	Z3 inverse(const Z3&);
	
	//IO operations. load() reduces, so a structure table written by BPtab
	//with honest Z_(p) integers in it arrives already reduced mod p. Note
	//this can turn a stored coefficient into 0, leaving an explicit zero
	//term in the vector -- BP_Op's reduce_*_mod_I helpers are what drop
	//those, and they must be run over anything loaded from a BPtab file.
	string output(Z3);
	void save(const Z3&, std::iostream&);
	Z3 load(std::iostream&);
	
	//the canonical representative: the residue mod p
	Z3 normalize(const Z3&);
};
