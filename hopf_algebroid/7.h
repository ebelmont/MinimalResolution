//embed into a cofree comodule
#include<omp.h>
template<typename ring, typename algebroid>
template<typename degree_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::embed2cofree(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, curtis_table<ring> *table, std::vector<int> *gens, matrix<ring> *new_map, std::vector<int> *basis_order, std::iostream &tablefile){
	std::cout << X->rank() << std::flush;
	std::vector<int> gns;
	if(gens==NULL) gens = &gns;
	
	//initialize the map to a zero comodule
	inj->set2zero(X->rank());
	table->clear();
	cofree_comodule<algebroid,degree_type> res;
//	std::cout << "KK" << std::flush;
	//the set of basis
	std::vector<int> bas;
	for(int i=0; i<X->rank(); ++i)
		bas.push_back(i);
	
	//sort the basis under degree
	auto i_deg = [X] (int m){ 
		return cofree_comodule<algebroid,degree_type>::underlyingDeg(X->degree(m)); };
	std::cout << "sorting..." << std::flush;
	std::stable_sort(bas.begin(),bas.end(),[i_deg](int m,int n) { return i_deg(m)<i_deg(n); });
	
	if(basis_order != NULL)
		*basis_order = bas;
	
//	int current_deg = 0;
	//throttled -- see matrices_mem/4.h::gaussian; on mr_st's non-modeled
	//resolvor path this loop runs X->rank() times, large at degree-270 scale
	for(int j=0;j<X->rank();++j){
		if(j % 500 == 0)
			std::cout << "\r" << j << "/" << X->rank() << std::flush;
		//look at the i-th row of the current map
		int i = bas[j];

		//load the i-th row
		vectors<matrix_index, ring> irow =  adjoint(X, *gens, res.position_of_gens, i, 0);
		//trace what is inside the current vector
		vectors<matrix_index, ring> sc = moduleOper->singleton(i);
		//simplify the i-th row using the current table, saving the homotopy in sc
		matrix_index pos = table->simplify_to_led(res.rank(),X->rank(),irow,sc);

		if(pos!=curtis_table<ring>::Boundary){
			//this means the i-th row maps to a nontrivial element in the cofree one
			table->insert(pos,i,irow,sc);
			//update inj
			inj->insert(i,irow);
		}
		else{
			//we need a new generator
			gens->push_back(i);

			//record the current rank
			int old_rank = res.rank();

			//compute the new summand of the cofree comodule
			cofree_comodule<algebroid,degree_type> new_sumd = adjoint(X,i,NULL,old_rank);

			//add the new summand
			res.direct_sum(new_sumd);

			//add new summand to the maps
			irow = adjoint(X, *gens, res.position_of_gens, i, 0);
			inj->insert(i,irow);

			//update the table
			table->insert(old_rank,i,irow,sc);
		}
	}
//	table->save(tablefile);
	return res;
}

//embed into a cofree one, using a model
template<typename ring, typename algebroid>
template<typename degree_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::embed2cofree_modeled(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, std::vector<int> *gens, matrix<ring>* new_map){
	//initialize
	std::cout << "a" << std::flush;
	inj->set2zero(X->rank());
	std::cout << "o" << std::flush;
	cofree_comodule<algebroid,degree_type> res;
	new_map->clear();
	cofree_comodule<algebroid,degree_type> new_sumd; 

	for(int i = 0; i<(int)gens->size(); ++i){
		//throttled -- see matrices_mem/4.h::gaussian
		if(i % 500 == 0)
			std::cout << "\r" << i << "/" << gens->size() << std::flush;
		int old_rank = res.rank();
		//get the new summand
		new_sumd = adjoint(X,(*gens)[i],new_map,0);
		//insert the new summand
		res.direct_sum(new_sumd);
		inj->direct_sum(*new_map,old_rank);
		new_map->clear();
	}
	return res;
}

//embed into a cofree one, using a model
template<typename ring, typename algebroid>
template<typename degree_type>
cofree_comodule<algebroid,degree_type> Hopf_Algebroid<ring,algebroid>::embed2cofree_modeled(const CoModule<algebroid,degree_type> *X, matrix<ring> *inj, std::vector<int> *gens){
	//initialize
	cofree_comodule<algebroid,degree_type> res;
	cofree_comodule<algebroid,degree_type> new_sumd;
	std::cout << "a" << gens->size() << std::flush;
	auto sum_timer = std::unique_ptr<DiagTimer>(new DiagTimer("summands"));
	//construct the new comodule
	//throttled -- see matrices_mem/4.h::gaussian; this is the
	//embed2cofree_modeled overload resolvor_modeled actually calls, so this
	//loop runs once per resolution stage
	for(int i = 0; i<(int)gens->size(); ++i){
		if(i % 500 == 0)
			std::cout << "\r" << i << "/" << gens->size() << "   " << std::flush;
		//get the new summand
		new_sumd = adjoint(X, (*gens)[i], NULL,0);
		//insert the new summand
		res.direct_sum(new_sumd);
	}

	//construct the new map -- likely the largest loop bound in the whole
	//resolve() pipeline (X->rank() iterations). Parallelized in two
	//stages: (1) gather X->coaction(i) for all i serially first, since
	//X->coaction can bottom out in a matrix_file-backed table with shared
	//stream state (e.g. BP_Op::delta -> delta_table) -- racing that would
	//serialize on its internal lock at best, corrupt at worst; see
	//hopf_algebroid/10.h's quotient comment for the reference hazard
	//analysis this mirrors. (2) run the pure per-row math
	//(adjoint_from_coaction) in parallel via matrix<ring>::construct_parallel,
	//with any conversion table it needs (e.g. BP_Op's R2L_table, also
	//possibly matrix_file-backed) materialized for the duration via
	//with_materialized_conversion_tables.
	sum_timer.reset();
	int n = X->rank();
	std::vector<vectors<matrix_index,algebroid>> coas(n);
	{
		DiagTimer _dt("gather_coaction");
		for(int i=0; i<n; ++i){
			if(i % 500 == 0)
				std::cout << "\r" << i << "/" << n << std::flush;
			coas[i] = X->coaction(i);
		}
	}

	std::vector<int> const& gvec = *gens;
	std::vector<uint32_t> const& pos = res.position_of_gens;
	std::function<vectors<matrix_index,ring>(int)> map_rows = [this,&gvec,&pos,&coas](int i){
		return this->template adjoint_from_coaction<degree_type>(gvec, pos, coas[i], 0);
	};
	{
		DiagTimer _dt("map_rows", std::max(1, omp_get_max_threads()));
		with_materialized_conversion_tables([&](){
			inj->construct_parallel(n, map_rows);
		});
	}
	
	return res;
}
