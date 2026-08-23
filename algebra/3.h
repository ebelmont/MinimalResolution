template<typename R>
class RingOp : virtual public AbGroupOp<R>{
public:
	//multliplication
	virtual R multiply(const R&, const R&)=0;
	  
	//the unit map from the integers
	virtual R unit(int)=0;
	  
	//if an element is invertible
	virtual bool invertible(const R&)=0;
	  
	//the inverse of an invertible element
	virtual R inverse(const R&)=0;

	//the valuation of an element with respect to the ring's distinguished
	//prime, for rings with such a notion (e.g. p-adic-like rings). Default:
	//every nonzero element is a unit, i.e. valuation 0 -- correct for plain
	//fields, overridden by rings that actually have non-unit elements.
	virtual unsigned valuation(R){ return 0; }

	//divide by (the distinguished prime)^n, for rings with a valuation.
	//Default is the identity, since a generic ring has no such structure
	//and this is only ever meaningfully called on a ring that overrides it.
	virtual R divide(const R& x, int){ return x; }

	//power of an elememnt
	R power(R x, unsigned n){
		if(n==0) return unit(1);
		if(n==1) return x;
		auto y = power(x,n/2);
		y = multiply(y,y);
		if(n%2 == 1) return multiply(y,x);
		else return y;
	}
};
