//Z3.h
#pragma once
#include"algebra.h"
#include"Fp.h"

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
	bool invertible(const Z3&);
	//inverse of an invertible element
	Z3 inverse(const Z3&);
};
