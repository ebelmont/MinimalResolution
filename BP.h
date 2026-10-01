//BP.h
#pragma once
#include"Z3.h"
#include"polynomial.h"
#include"matrices.h"
#include"mon_index.h"
#include"hopf_algebroid.h"
#include"BPQ.h"
#include<map>

//BP_* = Zp[v_1,v_2,\dots]
typedef polynomial<Z3> BP;
//BP_*BP = BP_*[t_1,t_2,\dots]
typedef polynomial<BP> BPBP;
typedef polynomial<BPBP> BPBPBP;
    
//operations for BP, with BPBP as a polynomial algebra over BP via the Right unit.
class BP_Op;

//operations on BP_*BP    
class BPBP_Op : public PolynomialOp_Para<BP>{
public:
	BPBP_Op(BP_Op*);
};
    
//operations on BPBPBP
class BPBPBP_Op : public PolynomialOp_Para<BPBP>{
public:
	BPBPBP_Op(BPBP_Op*);
};
    
//structur of Hopf algebroid on (BP_*,BP_*BP)
class BP_Op : virtual public Hopf_Algebroid<BP,BPBP>, public PolynomialOp_Para<Z3>{
public:
	//ring operations
	Z3_Op *Z3_oper;
	BPBP_Op BPBP_opers;
	BPBPBP_Op BPBPBP_opers;
        
	//module operations
	ModuleOp<matrix_index,BP> BPMod_opers;
	ModuleOp<matrix_index,BPBP> BPBPMod_opers;
	ModuleOp<matrix_index,Z3> Z3Mod_opers;

	//index of monomials
	monomial_index mon_index;
	
	//The height n of the invariant ideal I_n = (p, v_1, ..., v_{n-1}) that
	//this operator works modulo.
	//
	//  height == 0  -- the classical case. The base ring is BP_* itself and
	//                  every reduction below is the identity, so nothing
	//                  changes for mr_BP / BPtab / the existing comodules.
	//  height >= 1  -- the base ring is BP_*/I_n = F_p[v_n, v_{n+1}, ...] and
	//                  the algebroid is Gamma(n) = BP_*BP/I_n. I_n is an
	//                  invariant ideal, so this really is a Hopf algebroid,
	//                  and for a comodule M with I_n M = 0 the change-of-rings
	//                  isomorphism
	//                      Ext_{BP_*BP}(BP_*, M) = Ext_{Gamma(n)}(BP_*/I_n, M)
	//                  says resolving M over it computes the same thing. That
	//                  is what lets comodules like BP_*(S/p) = BP_*/p -- free
	//                  over BP_*/I_1 but NOT over BP_* -- go through the same
	//                  machinery. See docs/GENERAL_COMODULES.md.
	//
	//Two independent things make up "mod I_n", and only the second lives
	//here: the coefficients become F_p (that is Z3_mod_p_Op, handed to the
	//constructor), and the monomials v_1,...,v_{n-1} die (that is the
	//reduce_*_mod_I helpers below, applied to the structure tables as they
	//are loaded). Killing those monomials is reduction modulo a MONOMIAL
	//ideal, so the reduced elements are closed under + and *: once the
	//tables and the comodule's coaction are reduced, every product the
	//resolution engine forms is automatically reduced too, and no ring
	//operation needs overriding.
	int height;
	
	//set the height. Call this BEFORE initialize(), since the reduction is
	//applied to the structure tables while they are being read.
	void set_height(int n);
	
	//true if the v-monomial e involves one of v_1,...,v_{height-1}, i.e. if
	//it dies in BP_*/I_height. Always false at height 0.
	bool killed_v_monomial(exponent e);
	
	//Reduce mod I_height. Which helper to use depends on which of BPBP's two
	//exponent slots means what, and BPBP's layout is the opposite of what the
	//type name suggests -- see the warning at the top of comodules.cpp. In
	//the RIGHT-unit presentation (the one coaction entries and the etaL and
	//delta tables use) the outer exponent indexes the v_i via eta_R and the
	//inner BP's exponent indexes the t_i; in the LEFT-unit presentation (what
	//R2L produces and its table stores) it is the other way round.
	//
	//All four drop terms that become zero, which matters: a stray zero
	//coefficient would survive scalor_mult and mon_multiply (neither of them
	//filters) and could end up as the leading term handed to
	//PolyOp::inverse.
	//
	//a genuine BP_* element, exponent = v-monomial
	BP reduce_v_mod_I(const BP&);
	//a polynomial in the t_i, exponent = t-monomial: only the coefficients move
	BP reduce_t_mod_I(const BP&);
	//BPBP in the right-unit presentation: outer = v via eta_R, inner = t
	BPBP reduce_right_mod_I(const BPBP&);
	//BPBP in the left-unit presentation: outer = t, inner = v
	BPBP reduce_left_mod_I(const BPBP&);
        
private:
	//table for eta_L
	matrix<BP> *etaL_table;
        
    //table for delta
	matrix<BPBP> *delta_table;
	
	//table for change from right to left and switch the place of vi and ti
	matrix<BP> *R2L_table;
        
public:
	BP_Op(int maxdeg, Z3_Op*, matrix<BP> *etaL_mat, matrix<BPBP> *delta_mat, matrix<BP> *R2L_mat);
        
	//load the structure maps from files
	void load_etaL(string);
	void load_R2L(string);
	void load_delta(string);

	//left unit
	BPBP etaL(const BP&);
	
	//change the notation from right unit presentation to left unit presentation
	BPBP etaL(const BPBP&);
        
	//right unit
	BPBP etaR(const BP&);
        
	//change the notation from right to left
	BPBP R2L(const BPBP&);
        
	//divide by p^n
	BP divide_power_p(const BP&,int);
	BPBP divide_power_p(const BPBP&,int);
	//divide powers of v1
	BPBP divide_v1(const BPBP&,int);
        
	//change from polynomial notation to vector notation
	vectors<matrix_index, BP> algebroid2vector(const BPBP&,int);
	
	//change from polynomial notation to vector notation
	vectors<matrix_index, BP> algebroid2vector(BPBP&&,int);
    
	//change from vector to polynoimial notation
	BPBP vector2algebroid(const vectors<matrix_index, BP>&);
    
	//the co-multiplication
	vectors<matrix_index, BPBP> delta(matrix_index);

	//number of generators below a degree
	unsigned ranksBelowDeg(unsigned); 
        
	//initilaizaiotn
	void initialize(string etaL_filename, string R2L_filename, string delta_filename);
        
	//lift elements in F3 to BP_*
	BP lift(F3);
	
	//lift vectors over F3 to vectors over BP_*
	vectors<matrix_index,BP> lift(const vectors<matrix_index,F3>&);
	
	//make the structure tables
	void make_tables(int maxVar, string R2Lfilename, string deltafilename, string etaL_filename, string R2L_filename, string delta_filename, std::ostream &outputfile=std::cout);
	
	//return the element t_1, built directly rather than from the structure
	//tables: outer exponent 0 (no v's), inner exponent the first t-slot
	BPBP t1();
	//return the element h0
	BPBP h0();
	//return the element v1
	BP v1();
	//return the top thetas on the Moore spectrum
	std::vector<BPBP> thetas();
};
    
//cofree comodule for BP_*BP
typedef cofree_comodule<BPBP,int> FreeBPCoMod;
//generic comodules
typedef comodule_generic<BPBP,int> BPCoMod_generic;
    
