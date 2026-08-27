#pragma once
#include<cstdlib>

//per-site OpenMP ablation switch, controlled by the OMP_SITE_MASK env var
//(bit K enables site K; unset = all on, i.e. current behavior). Purely a
//diagnostic tool for measuring which of the parallel sites re-enabled after
//00405f4 actually help vs. hurt at real scale -- see NOTES.md 2026-08-27.
//Sites: 0=sum (modules/8.h), 1=maps_to_p (matrices/5.h),
//2=mon_multiply (polynomial/5.h), 3=construct_parallel (matrices/2.h),
//4=row_reduction (matrices_mem/4.h), 5=toDense (modules/2.h),
//6=add2Dense (modules/2.h), 7=scalor_mult (modules/3.h)
inline bool omp_site_enabled(int site){
	static long mask = [](){
		const char *e = std::getenv("OMP_SITE_MASK");
		return e ? std::strtol(e, nullptr, 0) : ~0L;
	}();
	return (mask >> site) & 1L;
}
