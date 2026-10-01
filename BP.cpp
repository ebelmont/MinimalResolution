//BP.cpp
#include"BP.h"
#include"mon_index.h"
#include<ios>

//constructor
BP_Op::BP_Op(int maxdeg, Z3_Op *Z3_op, matrix<BP> *etaL_mat, matrix<BPBP> *delta_mat, matrix<BP> *R2L_mat) :
ModuleOp<exponent,Z3>(Z3_op), PolynomialOp_Para<Z3>(Z3_op), BPBP_opers(this), BPBPBP_opers(&BPBP_opers), BPMod_opers(this), BPBPMod_opers(&BPBP_opers), Z3Mod_opers(Z3_op), mon_index(maxdeg){
	mon_index.init_mon_array();
	maxDeg = maxdeg;
	Hopf_Algebroid<BP,BPBP>::ringOper = this;
	algebroidRingOper = &BPBP_opers;
    
	moduleOper = &BPMod_opers;
	algebroidModuleOper = &BPBPMod_opers;
    
	Z3_oper = Z3_op;
	//the classical case unless set_height says otherwise
	height = 0;
	etaL_table = etaL_mat;
	delta_table = delta_mat;
	R2L_table = R2L_mat;
}

//constructor
BPBP_Op::BPBP_Op(BP_Op *BP_op) : ModuleOp<exponent,BP>(BP_op), PolynomialOp_Para<BP>(BP_op){}

//constrctor
BPBPBP_Op::BPBPBP_Op(BPBP_Op *BPBP_op) : ModuleOp<exponent,BPBP>(BPBP_op), PolynomialOp_Para<BPBP>(BPBP_op){}

//set the height
void BP_Op::set_height(int n){
	if(n < 0){
		std::cerr << "BP_Op::set_height: negative height " << n << ", using 0\n" << std::flush;
		n = 0;
	}
	height = n;
}

//true if the v-monomial e involves one of v_1,...,v_{height-1}
bool BP_Op::killed_v_monomial(exponent e){
	for(int i=1; i<height; ++i)
		if(xnVal(e,i) != 0) return true;
	return false;
}

//reduce a genuine BP_* element mod I_height
BP BP_Op::reduce_v_mod_I(const BP &x){
	if(height == 0) return x;
	BP result;
	for(auto &tm : x.dataArray){
		//monomials involving v_1,...,v_{height-1} lie in I_height
		if(killed_v_monomial(tm.ind)) continue;
		//and so do coefficients divisible by p, once the coefficient ring is F_p
		Z3 c = Z3_oper->normalize(tm.coeficient);
		if(Z3_oper->isZero(c)) continue;
		result.push({tm.ind, c});
	}
	return result;
}

//reduce a polynomial in the t_i: the t_i survive, only the coefficients move
BP BP_Op::reduce_t_mod_I(const BP &x){
	if(height == 0) return x;
	BP result;
	for(auto &tm : x.dataArray){
		Z3 c = Z3_oper->normalize(tm.coeficient);
		if(Z3_oper->isZero(c)) continue;
		result.push({tm.ind, c});
	}
	return result;
}

//reduce a BPBP written in the right-unit presentation: outer = v, inner = t
BPBP BP_Op::reduce_right_mod_I(const BPBP &x){
	if(height == 0) return x;
	BPBP result;
	for(auto &tm : x.dataArray){
		if(killed_v_monomial(tm.ind)) continue;
		BP c = reduce_t_mod_I(tm.coeficient);
		if(c.size() == 0) continue;
		result.push({tm.ind, c});
	}
	return result;
}

//reduce a BPBP written in the left-unit presentation: outer = t, inner = v
BPBP BP_Op::reduce_left_mod_I(const BPBP &x){
	if(height == 0) return x;
	BPBP result;
	for(auto &tm : x.dataArray){
		//the outer slot is a t-monomial and never dies
		BP c = reduce_v_mod_I(tm.coeficient);
		if(c.size() == 0) continue;
		result.push({tm.ind, c});
	}
	return result;
}

//Load a structure table, reducing each row mod I_height as it is read.
//
//The rows of these tables are indexed by mon_index position (row i is the
//value on mon_array[i]) and carry raw-exponent-indexed values; matrix_index
//and exponent are the same 32-bit type, which is why a row of a matrix<BP>
//is literally a BPBP. matrix::construct walks the rows in order and only
//calls clear()/set_rank()/insert() on the TARGET, so reading the source
//stream sequentially from inside row_rule is safe, and it works for the
//file-backed matrix_file backend, whose update_all is unimplemented.
//
//row_kills says whether row i itself dies: for the tables indexed by
//v-monomials a row whose monomial involves v_1,...,v_{height-1} is the value
//on something that is 0 in BP_*/I_height. Those rows are never consulted once
//every input is reduced, but zeroing them keeps the invariant "everything in
//this table lies in Gamma(height)" literally true.
template<typename R>
static void load_reduced(matrix<R> *table, std::iostream &reader, int rk,
                          ModuleOp<matrix_index,R> *modoper,
                          std::function<vectors<matrix_index,R>(const vectors<matrix_index,R>&)> reduce,
                          std::function<bool(int)> row_kills){
	std::function<vectors<matrix_index,R>(int)> row_rule = [&](int i){
		auto row = modoper->load(reader);
		if(row_kills(i)) return modoper->zero();
		return reduce(row);
	};
	table->construct(rk, row_rule);
}

//load etaL table
void BP_Op::load_etaL(string filename){  
	std::cout << "loading etaL table from " << filename << "...\n" << std::flush;
	std::fstream reader(filename, std::ios::in | std::ios::binary);
	if(height == 0)
		etaL_table->load(reader, mon_index.number_of_all_mons()); 
	else{
		//values are in the right-unit presentation; rows are indexed by v-monomials
		std::function<vectors<matrix_index,BP>(const vectors<matrix_index,BP>&)> reduce =
		    [this](const BPBP &x){ return reduce_right_mod_I(x); };
		std::function<bool(int)> row_kills = [this](int i){
			return killed_v_monomial(mon_index.mon_array[i]); };
		load_reduced<BP>(etaL_table, reader, mon_index.number_of_all_mons(),
		                  &BPMod_opers, reduce, row_kills);
	}
	std::cout << "etaL data loaded\n" << std::flush;
}

//load R2L talbe
void BP_Op::load_R2L(string filename){   
	std::cout << "loading R2L table...\n" << std::flush;
	std::fstream reader(filename, std::ios::in | std::ios::binary);
	if(height == 0)
		R2L_table->load(reader, mon_index.number_of_all_mons());
	else{
		//values are in the left-unit presentation; rows are indexed by v-monomials
		std::function<vectors<matrix_index,BP>(const vectors<matrix_index,BP>&)> reduce =
		    [this](const BPBP &x){ return reduce_left_mod_I(x); };
		std::function<bool(int)> row_kills = [this](int i){
			return killed_v_monomial(mon_index.mon_array[i]); };
		load_reduced<BP>(R2L_table, reader, mon_index.number_of_all_mons(),
		                  &BPMod_opers, reduce, row_kills);
	}
	std::cout << "R2L data loaded\n" << std::flush;
}

//load delta table
void BP_Op::load_delta(string filename){    
	std::cout << "loading delta table...\n" << std::flush;
	std::fstream reader(filename, std::ios::in | std::ios::binary);
	if(height == 0)
		delta_table->load(reader, mon_index.number_of_all_mons());
	else{
		//a row is a sum of (t-monomial index, BPBP in the right-unit
		//presentation) pairs; rows are indexed by t-monomials, so no row dies
		std::function<vectors<matrix_index,BPBP>(const vectors<matrix_index,BPBP>&)> reduce =
		    [this](const vectors<matrix_index,BPBP> &x){
			vectors<matrix_index,BPBP> result;
			for(auto &tm : x.dataArray){
				BPBP c = reduce_right_mod_I(tm.coeficient);
				if(BPBP_opers.isZero(c)) continue;
				result.push({tm.ind, c});
			}
			return result;
		};
		std::function<bool(int)> row_kills = [](int){ return false; };
		load_reduced<BPBP>(delta_table, reader, mon_index.number_of_all_mons(),
		                    &BPBPMod_opers, reduce, row_kills);
	}
	std::cout << "delta data loaded\n" << std::flush;
}

//the left unit
BPBP BP_Op::etaL(const BP &x){
	std::function<BPBP(const Z3&, const BPBP&)> scalor_mult = [this] (const Z3 &r, const BPBP &fm){
		return BPBP_opers.scalor_mult(this->monomial(0,r),fm); };
	auto xvec = mon_index.poly2vec(x, this);
	return etaL_table->maps_to(xvec, scalor_mult, &BPBP_opers);
}

//change the left notation to the right notation
BPBP BP_Op::etaL(const BPBP &x){
	std::function<BPBP(const BP&, const BPBP&)> scalor_mult = [this](const BP &r, const BPBP &fm){
		return BPBP_opers.scalor_mult(r,fm); };
	auto xvec = mon_index.poly2vec(x,&BPBP_opers);
	return etaL_table->maps_to(xvec,scalor_mult, &BPBP_opers);
}

//the right unit, vn is in the outer
BPBP BP_Op::etaR(const BP &x){
	std::function<BP(const Z3&)> tfm = [this](const Z3& r){
		return this->monomial(0,r); };
	return this->termwise_operation(tfm,x);
}

//change algebroid to a vector
vectors<matrix_index, BP> BP_Op::algebroid2vector(const BPBP& x, int shift){
	std::function<matrix_index(exponent)> rd = [this,shift](exponent e){
		return mon_index.mon_index[e]+shift; };
	return BPBP_opers.re_index(rd,std::move(R2L(x)));
}

//change right notation to the left notation and switch ti to the outer
BPBP BP_Op::R2L(const BPBP &x){
	std::function<BPBP(const BP&, const BPBP&)> scalor_mult = [this](const BP &r, const BPBP &fm){
		std::function<BP(const Z3&)> tfm = [this](const Z3& r0){
			return this->monomial(0,r0); };
		auto r1 = termwise_operation(tfm,r);
		return BPBP_opers.multiply(r1,fm); };
	auto xvec = mon_index.poly2vec(x,&BPBP_opers);
	return R2L_table->maps_to(xvec,scalor_mult, &BPBP_opers);
}

//change back
BPBP BP_Op::vector2algebroid(const vectors<matrix_index, BP>&x){
    std::cerr << "not implemented!" ;
    return BPBP_opers.zero();
}
    
//the coaction
vectors<matrix_index, BPBP> BP_Op::delta(matrix_index n){
    return delta_table->find(n);
}

//number of generators
unsigned BP_Op::ranksBelowDeg(unsigned n){ 
	return mon_index.ranksBelow[n]; }

//initialize
void BP_Op::initialize(string etaL_filename, string R2L_filename, string delta_filename){
	load_etaL(etaL_filename);
	load_R2L(R2L_filename);
	load_delta(delta_filename);
    
	//initialize the degree of generators
	std::function<int(matrix_index)>  cofree_degs = [this](matrix_index n){
		exponent e = mon_index.mon_array[n];
		return total_deg(e,xnDeg);
	};
        
	//initialize the data for a cofree coalgebra
	this->init_cofree_data(cofree_degs);
}

//make structure tables
void BP_Op::make_tables(int maxVar, string R2Lfilename, string deltafilename, string etaL_filename, string R2L_filename, string delta_filename, std::ostream &outputfile){
	//the values on the generators
	std::vector<BPBP> R2L_gen(maxVar+1);
	std::vector<BPBPBP> delta_gen(maxVar+1);
	std::vector<BPBP> etaL_gen(maxVar+1);
	
	//load the values on the generators. 
	//for etaR, At this point we use the left unit expression, and the vi are on the outer, ti on the innner
	std::fstream R2L_file(R2Lfilename, std::ios::in | std::ios::binary);
	R2L_gen[0] = BPBP_opers.unit(1);
	for(int i=1; i<=maxVar; ++i)
		R2L_gen[i] = BPBP_opers.load(R2L_file);
	R2L_file.close();
	std::cout << "etaR data loaded\n";
	for(int i=1; i<=maxVar; ++i)
		outputfile << "etaR(v" << i << ") = " << BPBP_opers.output(R2L_gen[i]) << "\n";
	
	//for delta, at this point, we use left unit expression, the vi in the middle, right ti in the outer, the left ti in the inner
	std::fstream delta_file(deltafilename, std::ios::in | std::ios::binary);
	delta_gen[0] = BPBPBP_opers.unit(1);
	for(int i=1; i<=maxVar; ++i)
		delta_gen[i] = BPBPBP_opers.load(delta_file);
	delta_file.close();
	std::cout << "delta data loaded\n";
	for(int i=1; i<=maxVar; ++i)
		outputfile << "delta(t" << i << ") = " << BPBPBP_opers.output(delta_gen[i]) << "\n";
	
	//compute the generators for etaL, in right unit expression,
	etaL_gen[0] = BPBP_opers.unit(1);
	for(int i=1; i<=maxVar; ++i){
		PolynomialOp_Para<BPBP> PBPBPoper(&BPBP_opers);
		polynomial<BPBP> etaL_rule = PBPBPoper.constant(BPBP_opers.monomial(vars(i)));
		//etR = -etaR(vn) + etaL(vn) in left expression , with vi in the outer
		auto etR = R2L_gen[i];
		etR = BPBP_opers.minus(etR);
		etR = BPBP_opers.add(BPBP_opers.monomial(vars(i)), etR);
		
		//lift etR to a polynomial
		std::function<BPBP(const BP&)> ct = [this](const BP& x){ 
			return BPBP_opers.constant(x); };
		auto tR = BPBP_opers.termwise_operation(ct,etR);
		
		//etaL(v_n) = -etaR(vn) + etaL(vn) + etaR(vn)
		etaL_rule = PBPBPoper.add(etaL_rule, tR);
		
		etaL_gen[i] = substitute(etaL_gen, etaL_rule, &BPBP_opers);
	}
	std::cout << "etaL data computed\n";
	for(int i=1; i<=maxVar; ++i)
		outputfile << "etaL(v" << i << ") = " << BPBP_opers.output(etaL_gen[i]) << "\n";
	
	//compute the etaL table
	std::function<BPBP(int)> etaL_gens = [&etaL_gen](int i){
		return etaL_gen[i]; };
	std::fstream etaL_tablefile(etaL_filename, std::ios::out | std::ios::binary);
	mon_index.substitution_table(etaL_gens, etaL_tablefile, &BPBP_opers);
	etaL_tablefile.close();
	
	//compute the R2L table
	std::function<BPBP(int)> R2L_gens = [&R2L_gen, this](int i){
		//switch vi to the outer
		std::function<BPBP(BPBP&&,BPBP&&)> merger = [this](BPBP &&x, BPBP &&y){
			return BPBP_opers.add(std::move(x),std::move(y)); };
		return swapping(R2L_gen[i], merger); };
	std::fstream R2L_tablefile(R2L_filename, std::ios::out | std::ios::binary);
	mon_index.substitution_table(R2L_gens, R2L_tablefile, &BPBP_opers);
	R2L_tablefile.close();
	
	//load the etaL table
	load_etaL(etaL_filename);
	//compute the delta table
	std::function<vectors<matrix_index, BPBP>(const BPBPBP&)> rule = [this](const BPBPBP& srf){
		//change from left unit notation to right unit notation
		std::function< std::pair<matrix_index,BPBP>(exponent,const BPBP&)> cf = [this](exponent e, const BPBP &w){
			return std::make_pair(mon_index.mon_index[e],std::move(etaL(w))); };
		 return BPBPBP_opers.termwise_operation(cf,srf); };
	std::function<void(const vectors<matrix_index,BPBP>&, std::iostream&)> outputer = [this](const vectors<matrix_index,BPBP>& x, std::iostream& writer){
		BPBPMod_opers.save(x, writer); };
	std::function<BPBPBP(int)> delta_gens = [&delta_gen](int i){
		return delta_gen[i]; };
	std::fstream delta_tablefile(delta_filename, std::ios::out | std::ios::binary);
	mon_index.substitution_table(delta_gens, delta_tablefile, &BPBPBP_opers, rule, outputer);
	delta_tablefile.close();
}

//lift elements in F3 to BP
BP BP_Op::lift(F3 x){
	return monomial(0,Z3_oper->lift(x));
}

//lift vectors over F3 to vectors over BP
vectors<matrix_index,BP> BP_Op::lift(const vectors<matrix_index,Fp>&v){
	vectors<matrix_index, BP> result;
	for(auto tm : v.dataArray)
		result.push({tm.ind, lift(tm.coeficient)});
	return result;
}

//the degree of comodules
template<>
int FreeBPCoMod::underlyingDeg(int i){
	return i; }
	
//add degrees
template<>
int FreeBPCoMod::add_degree(int a, int b){
	return a+b; }

//the element v1
BP BP_Op::v1(){
	BP v1 = this->singleton(1);
	return v1;
}
	
//the element t_1, built directly from the exponent encoding rather than from
//the structure tables. Outer exponent 0 (no v's), inner exponent the first
//t-slot -- see the warning at the top of comodules.cpp about which slot is
//which. This is bitwise equal to h0() at height 0.
BPBP BP_Op::t1(){
	BP inner = monomial(singleVar(1,1), Z3_oper->unit(1));
	return BPBP_opers.monomial(0, inner);
}

//return the element h0 = (etaR(v1) - etaL(v1))/p
BPBP BP_Op::h0(){
	if(mon_index.max_degree<=1){
		std::cerr << "out of range for h0";
		return BPBP_opers.zero();
	}
	//In characteristic p the defining formula computes 0/p: eta_L(v_1) and
	//eta_R(v_1) differ by p*t_1, so their difference vanishes before
	//divide_power_p ever sees it, and dividing by p is not an operation
	//BP_*/I_n has. The answer is t_1 either way, so build it directly.
	if(height > 0){
		BPBP h0 = t1();
		std::cout << "h0=" << BPBP_opers.output(h0) << "\n";
		return h0;
	}
	auto dv1 = BPBP_opers.add(BPBP_opers.minus(etaL(v1())), etaR(v1()));
	BPBP h0 = divide_power_p(dv1,1);
	std::cout << "h0=" << BPBP_opers.output(h0) << "\n";
	return h0;
}

//return the top thetas on the Moore spectrum
std::vector<BPBP> BP_Op::thetas(){
	std::vector<BPBP> theta(10);
	if(mon_index.max_degree<=12){
		std::cerr << "out of range for beta1";
		return theta;
	}
	//v2
	BP v2 = monomial(vars(2));
	std::cout << "v2=" << output(v2) << "\n";
	//d(v2)
	BPBP dv2 = BPBP_opers.add(BPBP_opers.minus(etaL(v2)), etaR(v2));
	//beta1 = d(v2)/v1
	BPBP beta1 = divide_v1(dv2,1);
	std::cout << "beta1=" << BPBP_opers.output(beta1) << "\n";
	theta[2] = beta1;

	if(mon_index.max_degree<=28){
		std::cerr << "out of range for beta2";
		return theta;
	}
	//v2^2
	BP v22 = multiply(v2,v2);
	//d(v2^2)
	BPBP dv22 = BPBP_opers.add(BPBP_opers.minus(etaL(v22)), etaR(v22));
	//beta2 = d(v2^2)/v1
	BPBP beta2 = divide_v1(dv22,1);
	std::cout << "beta2=" << BPBP_opers.output(beta2) << "\n";
	theta[3] = beta2;

	if(mon_index.max_degree<=36){
		std::cerr << "out of range for beta3/3";
		return theta;
	}
	//v2^3
	BP v23 = multiply(v22,v2);
	//d(v2^3)
	BPBP dv23 = BPBP_opers.add(BPBP_opers.minus(etaL(v23)), etaR(v23));
	//beta3/3 = d(v2^3)/v1^3
	BPBP beta33 = divide_v1(dv23,3);
	std::cout << "beta3/3=" << BPBP_opers.output(beta33) << "\n";
	theta[4] = beta33;


	if(mon_index.max_degree<=60){
		std::cerr << "out of range for beta4";
		return theta;
	}
	//v2^4
	BP v24 = multiply(v23,v2);
	//d(v2^4)
	BPBP dv24 = BPBP_opers.add(BPBP_opers.minus(etaL(v24)), etaR(v24));
	//beta4 = d(v2^4)/v1
	BPBP beta4 = divide_v1(dv24,1);
	std::cout << "beta4=" << BPBP_opers.output(beta4) << "\n";
	theta[5] = beta4;


	if(mon_index.max_degree<=76){
		std::cerr << "out of range for beta5";
		return theta;
	}
	//v2^5
	BP v25 = multiply(v24,v2);
	//d(v2^5)
	BPBP dv25 = BPBP_opers.add(BPBP_opers.minus(etaL(v25)), etaR(v25));
	//beta5 = d(v2^5)/v1
	BPBP beta5 = divide_v1(dv25,1);
	std::cout << "beta5=" << BPBP_opers.output(beta5) << "\n";
	theta[6] = beta5;

	if(mon_index.max_degree<=84){
		std::cerr << "out of range for beta6/3";
		return theta;
	}
	//v2^6
	BP v26 = multiply(v25,v2);
	//d(v2^6)
	BPBP dv26 = BPBP_opers.add(BPBP_opers.minus(etaL(v26)), etaR(v26));
	//beta6/3 = d(v2^6)/v1^3
	BPBP beta63 = divide_v1(dv26,3);
	std::cout << "beta6/3=" << BPBP_opers.output(beta63) << "\n";
	theta[7] = beta63;




	return theta;
}

//divide by p^n
BP BP_Op::divide_power_p(const BP& x, int n){
	std::function<Z3(const Z3&)> rl = [this,n](const Z3 &r){
		return Z3_oper->divide(r,n); };
	return this->termwise_operation(rl,x);
}

//divide by p^n
BPBP BP_Op::divide_power_p(const BPBP& a, int n) {
	std::function<BP(const BP&)> rl = [this,n](const BP& r){ 
		return divide_power_p(r,n);};
	return BPBP_opers.termwise_operation(rl,a);
}

//divide by v1^n
BPBP BP_Op::divide_v1(const BPBP& x,int n){
	//mod 3 reduction
	BPBP x1;
	for(auto tm : x.dataArray){
		BP a1;
		for(auto am : tm.coeficient.dataArray) {
			if(am.coeficient%3==1){
				a1.push({am.ind, 1});
			}
			else if(am.coeficient%3==2){
				a1.push({am.ind, 2});
			}
		}
		if(!isZero(a1))
			x1.push({tm.ind,a1});
	}
	
	//divide by v1^n
	std::function<exponent(exponent)> rl = [this,n](exponent e) { 
		auto nt = unpack(e);
		if(nt[0]<n) std::cerr << "not v1-divisible";
		nt[0] -= n;
		return pack(nt.data());
	};
	return BPBP_opers.re_index(rl,x1);
}
