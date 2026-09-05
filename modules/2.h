//data structure for vectors with entries in a ring R
template<typename index, typename R>
class vectors{
public:
	//the terms in a vector, including an index and a coeficient
	class term { 
	public:
		//the index of the term
		index ind;
		//the coeficient of the term
		R coeficient; 
		
		//the default constructor
		term(){}
		
		//constructor with given index and coeficient
		term(index i, const R &c) {
			ind = i; coeficient = c; }
			
		//constructor which is destructive for the inputs
		term(index i, R &&c) {
			ind=i; coeficient=std::move(c); }
				
		//comparison operator compares the index of the terms
		bool operator<(const term &y) { return ind<y.ind; }	
	};
	
	//default constructor
	vectors(){}
	
	//the data for the vector, which is an array of terms
	std::vector<term> dataArray; 
	
	//clear the contents in a vector
	inline void clear(){
		dataArray.clear(); }
		
	//the size of the vector
	inline unsigned size() const{
		return dataArray.size();}
			
	//add a new term in the back
	inline void push(const term &tm){
		dataArray.push_back(tm); } 
				
	//add a new term, destructive in the argument
	inline void push(term &&tm){
		dataArray.push_back(std::move(tm)); }
					
	//sort the terms in an unordered array
	virtual void sort(){
		std::sort(dataArray.begin(), dataArray.end()); }
						
	//move the data from an array of terms
	void get_data(vectors<index,R> &&y){ 
		dataArray = std::move(y.dataArray); }
							
	//direct sum of two vectors
	virtual void direct_sum(const vectors<index,R> &y, unsigned rank){
		for(auto &tm:y.dataArray){
			term newTerm(tm.ind + rank, tm.coeficient);
			push(std::move(newTerm));
		}
	}
							
	//the inverse operation of direct sum: split a vector into two vectors, with given splitting point
	vectors<index,R> un_direct_sum(unsigned rank){
		unsigned i;
		//find the term at which to split
		for(i=0; i<size(); ++i)
			if(dataArray[i].ind>=rank) break;
									
		//construct the new vector with the tail
		vectors<index,R> result;
		for(unsigned k=i; k<size(); ++k){
			term newTerm(dataArray[k].ind - rank, std::move(dataArray[k].coeficient));
			result.push(std::move(newTerm));
		}
		//truncate the tail from the old vector
		dataArray.resize(i);
		return result;
	}
	
	//change into dense vector, 0 must transform to zeo
	std::vector<R> toDense(unsigned rank){
		std::vector<R> result(rank);
		if(size()==0) return result;
		if(dataArray[size()-1].ind >= rank)
			std::cerr << "tout:" << dataArray[size()-1].ind << "/" << rank << std::flush;
		for(unsigned i=0; i<size(); ++i)
			result[(unsigned) dataArray[i].ind] = dataArray[i].coeficient;
		return result;
	}
	
	//change back from dense vector
	void deDense(std::vector<R> const& y, RingOp<R> *ringop, unsigned start){
		clear();
		for(unsigned i = start; i<y.size(); ++i)
			if(!ringop->isZero(y[i])){
				term newTerm((index) i, y[i]);
				push(newTerm);
			}
	}
	
	//add to a dense vector
	void add2Dense(R c, std::vector<R> &y, RingOp<R> *ringop){
		if(size()==0) return;
		//mirrors scalor_mult's own early-out (modules/3.h): when the scalar
		//is zero, every one of this row's terms would multiply out to zero
		//and then no-op-add back into `y` -- but without this check, that
		//"no-op" still executes a real ringop->add(y[ind], zero) per term,
		//and add()'s zero-branch returns its non-zero argument BY VALUE
		//from a const-ref, forcing a full copy-construction of whatever
		//(possibly multi-term) ring element already sits at y[ind] -- for
		//every term, on every call, useful or not. Measured directly (see
		//NOTES.md's gaussian_workvisits entry): the vast majority of calls
		//have c==0 (rows with no actual interaction with this pivot), so
		//skipping the whole loop here removes that vast majority of the
		//accumulate phase's cost in one O(1) check, same as the old
		//take_away/scalor_mult path already did.
		if(ringop->isZero(c)) return;
		if(dataArray[size()-1].ind >= y.size())
			std::cerr << "aout:" << dataArray[size()-1].ind << "/" << y.size() << std::flush;
		//NOTE: tried `std::move`-ing both arguments here to select `add`'s
		//destructive overload instead of the copying one (hypothesis: avoid
		//a needless copy-construction on every call) -- measured a small
		//net REGRESSION at degree 270 (see NOTES.md), not the expected win.
		//Root cause understood, not just abandoned: the per-*term* payload
		//at this level (a `term<exponent,Z3>` inside a BP polynomial) is a
		//trivial POD, so copy vs. move of a term is free either way; the
		//real remaining cost is `add()`'s own `result.dataArray.reserve(...)`
		//-- ONE fresh heap allocation per call for the *whole merged
		//vector's backing storage*, paid identically by both overloads,
		//not something overload choice touches. Actually eliminating it
		//would mean reusing a scratch buffer *inside* `add()` itself
		//(shared, generic code used everywhere in the codebase -- exactly
		//the "higher risk, touches shared code" fix flagged as a last
		//resort) -- not attempted this session, left for a future one that
		//wants to take that risk deliberately.
		for(unsigned i=0; i<size(); ++i){
			R nf = ringop->multiply(c, dataArray[i].coeficient);
			y[(unsigned) dataArray[i].ind] = ringop->add(y[(unsigned) dataArray[i].ind], nf);
		}
	}

	//like add2Dense, but for a nonzero scalar also appends every touched
	//column index (this row's own dataArray indices) into `touched` --
	//same single pass, no second scan. Lets a caller build a per-*row*
	//touched-column set sized to the actual useful work done on that row
	//(sum of the densities of only the pivots that had a nonzero
	//coefficient against it), instead of a whole block's shared union
	//that every row would otherwise pay for regardless of how much of it
	//actually applied to that particular row (see matrices_mem/4.h's
	//gaussian(), the per-row trailing loop). Zero-scalar calls still
	//short-circuit exactly as in add2Dense -- they touch nothing and add
	//nothing to `touched`.
	void add2DenseTracked(R c, std::vector<R> &y, RingOp<R> *ringop, std::vector<index> &touched){
		if(size()==0) return;
		if(ringop->isZero(c)) return;
		for(unsigned i=0; i<size(); ++i){
			R nf = ringop->multiply(c, dataArray[i].coeficient);
			y[(unsigned) dataArray[i].ind] = ringop->add(y[(unsigned) dataArray[i].ind], nf);
			touched.push_back(dataArray[i].ind);
		}
	}

	//change back from a dense vector into this, like deDense, but (a) moves
	//each nonzero entry's payload out of `y` instead of copying it, and (b)
	//resets exactly the entries it takes back to a default-constructed
	//(zero) R, combining deDense's job with the reset a caller reusing `y`
	//as a scratch buffer for the next vector would otherwise need as a
	//separate pass. Merges two already-sorted candidate lists of indices
	//that could be nonzero (rather than requiring the caller to pre-merge
	//them, or scanning the whole buffer) -- needed because `y`'s rank can
	//be much larger than either list (a real case in practice --
	//del_and_gaussian's column space can be several times the row count,
	//see `matrices_mem/4.h`'s `gaussian()`). Correct with duplicates within
	//or across the two lists: the first occurrence of an index that's
	//still nonzero extracts it and zeroes `y` there, so a later duplicate
	//sees zero and is skipped. Result terms come out in sorted order
	//(required: `a` and `b` must each be sorted on input).
	void deDenseFromSortedUnion(std::vector<R> &y, RingOp<R> *ringop, std::vector<index> const &a, std::vector<index> const &b){
		clear();
		size_t i = 0, j = 0;
		while(i < a.size() || j < b.size()){
			index idx;
			if(i < a.size() && (j >= b.size() || a[i] <= b[j])) idx = a[i];
			else idx = b[j];
			if(i < a.size() && a[i] == idx) ++i;
			if(j < b.size() && b[j] == idx) ++j;
			if(!ringop->isZero(y[idx])){
				term newTerm(idx, std::move(y[idx]));
				push(std::move(newTerm));
				y[idx] = R();
			}
		}
	}
};
