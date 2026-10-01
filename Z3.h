//Z3.h
#pragma once
#include"algebra.h"
#include"Fp.h"
#include <cstdint>
#include <vector>

//we use 64-bit unsigned integer to denote 3-adic numbers
typedef uint64_t Z3;

//the residue field
typedef Fp F3;

class Z3_Op : public virtual RingOp<Z3>{
public:
	//highest valuation Z/3^40 can distinguish from a genuine zero (MAX3 =
	//3^40, Z3.cpp): a value divisible by 3^40 is isZero()==true, so 39 is
	//the last valuation this ring can represent honestly. Deliberately NOT
	//raised by widening Z3 -- **the target of this computation is
	//Ext(BP_*/3^40)**, so Z/3^40 is the intended model, not a limitation to
	//engineer around (see NOTES.md "2026-09-10: correction"). NOTE: this
	//comment used to say "Ext(BP_*) agrees with Ext(BP_*/3^40) outside stem
	//0" -- that is FALSE and has misled at least one session into treating
	//correct output as a bug. The LES of 0 -> BP_* -3^40-> BP_* ->
	//BP_*/3^40 -> 0 gives Ext^s(BP_*/3^40) = Ext^s(BP_*)/3^40 (+)
	//Ext^{s+1}(BP_*), i.e. a whole extra copy of Ext^{s+1} at v0-valuation
	//40-k for a class of order 3^k. The resulting chart is roughly twice
	//the size of Ext(BP_*) and that is the correct answer. Callers that drive a
	//v0-power search/cancellation depth (pric, in BP_init.cpp) must clamp
	//to this constant instead of using resolution_length directly, or the
	//search speculates past what this ring can correctly cancel and
	//manufactures bogus surviving generators (NOTES.md "the v0>=40 bug").
	static constexpr unsigned MAX_VALUATION = 39;

	//the constructor
	Z3_Op();

	//the charateristic of the residue field
	int prime();
	
	//the operations on F3
	Fp_Op F3_opers;
	
	//the valuation
	unsigned valuation(Z3) override;

	//n-th power of p
	Z3 power_p(int n);

	//incremental cache for power_p: index n holds p^n once computed. get_tag/
	//naming's v0-power computation used to recompute power_p(j) from scratch
	//via naive recursion for every j, O(pric) redundant work per call --
	//power_p() below fills this incrementally (one multiply() per new
	//highest power ever requested) instead.
	std::vector<Z3> power_p_cache;

	// divide by p^n
	Z3 divide(const Z3&, int n) override;
	
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
