//matrix stored in memeries
template<typename ring>
class matrix_mem : public matrix<ring>{
	//the data for a matrix, consisting of an array of vectors
	std::vector<vectors<matrix_index,ring>> data;
public:
	//default constructor
	matrix_mem(){
		this->rank = 0; }
	
	matrix_mem(string){
		this->rank = 0; }
	
	//clear the contents
	void clear(){
		std::vector<vectors<matrix_index,ring>> emd;
		data.swap(emd);
		data.clear();
		this->rank = 0;
	}
	
	//find the n-th row
	vectors<matrix_index,ring> find(matrix_index n) const{
		return data[n];
	}
	
	//insert a new row
	void insert(matrix_index i, vectors<matrix_index,ring> const& x){
		if(i>=data.size()){
			data.resize(i+1);
			this->rank = i+1;
		}
		data[i] = x;
	}
	
	//set the rank
	void set_rank(unsigned n){
		this->rank = n;
		data.resize(n);
	}

	//construct with parallel algorithm. Safe here because set_rank()
	//pre-sizes `data` before the parallel loop, so each thread writes a
	//distinct pre-existing element -- no reallocation race (unlike
	//insert(), whose out-of-range resize is only safe single-threaded).
	//The base class declares this virtual (matrices/2.h) but had no
	//override anywhere -- every prior caller (e.g. quotient_p) was
	//silently running the serial fallback. Caller is responsible for
	//`rows` being safe to call concurrently for different indices -- see
	//hopf_algebroid/10.h's quotient() comment and
	//hopf_algebroid/7.h's embed2cofree_modeled for the reference hazard
	//analysis and how callers gather any unsafe reads serially first.
	void construct_parallel(int rk, std::function<vectors<matrix_index,ring>(int)> rows) override {
		clear();
		set_rank(rk);
		#pragma omp parallel for schedule(dynamic)
		for(int i=0; i<(int)this->rank; ++i)
			data[i] = rows(i);
	}

	//update all rows
	void update_all(std::function<void(vectors<matrix_index,ring>&,matrix_index)> action){
		for(unsigned i=0; i<data.size(); ++i)
			action(data[i],i);
	}
	
	//row transfermations, subtract the given row to cancel the given column	
	void row_reduction(matrix_index row, matrix_index col);
	
	//Gaussian ellimination, with given rows and colums
	void gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols);
	
	//delete some columns and then do Gaussion
	void del_and_gaussian(std::vector<std::pair<matrix_index,matrix_index>> const &row_cols, std::set<int> const &to_del){
		{ DiagTimer _dt("del_cols"); this->del_cols(to_del); }
		DiagTimer _dt("gaussian");
		gaussian(row_cols);
	}
	
	//unify a vector to make a given entry 1
	vectors<matrix_index,ring> unify(matrix_index n, vectors<matrix_index,ring> const &x);
	
	//ellimination a given entry
	void take_away(matrix_index n, vectors<matrix_index,ring> &x, vectors<matrix_index,ring> const &y);
};
