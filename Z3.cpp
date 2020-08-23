//Z3.cpp
#include"Z3.h"
#include<sstream>

constexpr uint64_t MAX3 = 12157665459056928801U; // 3^40 (largest power of 3 < 2^64)
constexpr uint64_t SQRT_MAX3 = 3486784401U; // 3^20
constexpr uint64_t MAX_INT64 = 18446744073709551615U; // 2^64-1

//the constructor
Z3_Op::Z3_Op() : F3_opers(3) {}

//the residual characteristic
inline int Z3_Op::prime(){ return 3; }

//addition
inline Z3 Z3_Op::add(const Z3& a, const Z3& b){
	uint64_t x = a % MAX3;
	uint64_t y = b % MAX3;
	if (y <= MAX_INT64 - x) { // no integer overflow issues, just use %
		return (x+y) % MAX3;
	}
	else if (x + y == 0) { // x + y = 2^64 (in Z)
		return (MAX_INT64 - MAX3) + 1; // MAX_INT64 = 2^64-1
	}
	else { // x + y in C++ is x + y - 2^64 in Z. We want x + y - 3^40 (in Z).
		return (x + y) + ((MAX_INT64 - MAX3) + 1);
	}
}

//addition
inline Z3 Z3_Op::add(Z3&& a, Z3&& b){
	uint64_t x = a % MAX3;
	uint64_t y = b % MAX3;
	if (y <= MAX_INT64 - x) { // no integer overflow issues, just use %
		return (x+y) % MAX3;
	}
	else if (x + y == 0) { // x + y = 2^64 (in Z)
		return (MAX_INT64 - MAX3) + 1; // MAX_INT64 = 2^64-1
	}
	else { // x + y in C++ is x + y - 2^64 in Z. We want x + y - 3^40 (in Z).
		return (x + y) + ((MAX_INT64 - MAX3) + 1);
	}
}

//multiplication
inline Z3 Z3_Op::multiply(const Z3& x, const Z3& y){
	// Write x = x0 + 3^20 x1, y = y0 + 3^20 y1. Then
	// x y = x0 y0 + 3^20 x0 y1 + 3^20 x1 y0 (mod 3^40), and x0 y0 < 3^40 < 2^64.
	uint64_t x0 = x % SQRT_MAX3;
	uint64_t x1 = x/SQRT_MAX3;
	uint64_t y0 = y % SQRT_MAX3;
	uint64_t y1 = y/SQRT_MAX3;
	// To calculate 3^20 x0 y1 mod 3^40, just need to know x0, y1 mod 3^20.
	uint64_t x0y1 = x0 * (y1 % SQRT_MAX3); // the product is < 3^40, so no overflow
	uint64_t y0x1 = y0 * (x1 % SQRT_MAX3); // the product is < 3^40, so no overflow
	return add(add(x0*y0, SQRT_MAX3 * (x0y1 % SQRT_MAX3)), SQRT_MAX3 * (y0x1 % SQRT_MAX3));
}

//unit
//EB: Is the only point of this to convert from int to Z3?
inline Z3 Z3_Op::unit(int x){ return x % MAX3; }

//check if zero
inline bool Z3_Op::isZero(const Z3& x){ return (x % MAX3) ==(Z3)0; }

//the zero element
inline Z3 Z3_Op::zero(){ return 0; }

//negation
inline Z3 Z3_Op::minus(const Z3& x){ return MAX3 - (x % MAX3); }

//output using hex form
string Z3_Op::output(Z3 x){
	std::stringstream res;
	res << std::hex << x;
	return res.str();
}

//write to a 8-byte. This would depend on the cpu archetecture
inline void Z3_Op::save(const Z3& x, std::iostream& writer){ writer.write((char*)&x, 8); } 

//load from a 8-byte
inline Z3 Z3_Op::load(std::iostream& reader){ 
	Z3 result;
	reader.read((char*)&result,8);
	return result;
}

//if invertible
bool Z3_Op::invertible(const Z3& x){ return x % prime() != 0; }

//the inverse, using a recursion formula
//EB: This appears to be equivalent to the first 63 terms of the Taylor series
//for 1/(1+y) with y = x-1
//An invertible element in Zp is a unit mod p. If p = 2, that means it's 1 + 2a_1 +...
//and this Taylor series centered at 1 converges 2-adically.
Z3 Z3_Op::inverse(const Z3& x){
	if (x % 3 == 1) {
		Z3 u = 1;
		for(int i=0;i<6;++i)
			u=multiply(u, add(2, minus(multiply(x,u))));
		return u;
	}
	else if (x % 3 == 2) { // Use Taylor series for 1/(1+y) with y = -x-1; need to negate
		Z3 u = 1;
		for(int i=0;i<6;++i)
			u=multiply(u, add(2, multiply(x,u)));
		return minus(u);
	}
	else {
		std::cout << "tried to invert non-invertible element " << x << "\n" << std::flush;
		throw "";

	}
}

//the 3-adic valuation
unsigned Z3_Op::valuation(Z3 x){
	//convention for valuation of 0
	//TODO: do we need to change max_val? I have no idea what's special about this number.
	static constexpr int max_val = 35536;
	if(isZero(x)) return max_val;
	
	//valuation of invertible is 0
	if(invertible(x)) return 0;
	//for non-invertible elements
	return valuation(x/prime()) + 1;
}

//p^n
Z3 Z3_Op::power_p(int n){
	if(n==0) return unit(1);
	Z3 p = unit(prime());
	if(n>0)
		return power_p(n-1)*p;
	//undefined if n<0
	std::cerr << "not integral!";
	return 0;
}

//divide by p^n. The missing term casued by truncation is guessed by the sign
//EB: I'm concerned because I have no idea what that comment ^ refers to
//FIXME
Z3 Z3_Op::divide(const Z3& x, int n){
	int64_t y=x;
	return y/(uint64_t)power_p(n);
}

//lift from residue field
Z3 Z3_Op::lift(F3 x){
    return x;
}
