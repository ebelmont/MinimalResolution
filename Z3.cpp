//Z3.cpp
#include"Z3.h"
#include<sstream>
#include<cmath>
#include<cassert>

//the constructor
Z3_Op::Z3_Op() : F3_opers(3) {}

//the residual characteristic
inline int Z3_Op::prime(){ return 3; }

//addition
inline Z3 Z3_Op::add(const Z3& x, const Z3& y){
	assert(x < pow(2,63) && y < pow(2,63)); // if this is hit, need to fix Z/2^64 issue
	return x+y;
}

//addition
inline Z3 Z3_Op::add(Z3&& x, Z3&& y){ return x+y; }

//multiplication
inline Z3 Z3_Op::multiply(const Z3& x, const Z3& y){
	assert(x < pow(2,32) && y < pow(2,32)); // if this is hit, need to fix Z/2^64 issue
	return x*y;
}

//unit
inline Z3 Z3_Op::unit(int x){ return x; }

//check if zero
inline bool Z3_Op::isZero(const Z3& x){ return x==(Z3)0; }

//the zero element
inline Z3 Z3_Op::zero(){ return 0; }

//negation
inline Z3 Z3_Op::minus(const Z3& x){ return -x;}

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
bool Z3_Op::invertible(const Z3& x){ return (x&(Z3)1) != 0; }

//the inverse, using a recursion formula
//EB: This appears to be equivalent to the Taylor series for 1/(1+y) with y = x-1
//TODO: Do we need to change the bound of 6?
Z3 Z3_Op::inverse(const Z3& x){
	Z3 u = 1;
	for(int i=0;i<6;++i)
		u=u*(2-x*u);
	return u;
}

//the 2-adic valuation
unsigned Z3_Op::valuation(Z3 x){
	//convention for 0
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

//divide by p. The missing term casued by truncation is guessed by the sign
Z3 Z3_Op::divide(const Z3& x, int n){
	int64_t y=x;
	return y/(int64_t)power_p(n);
}

//lift from residue field
Z3 Z3_Op::lift(F3 x){
    return x;
}
