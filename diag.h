//diag.h -- opt-in wall/cpu timing, off by default (DIAG_TIMING env var).
//
//A DiagTimer times its own scope. Timers opened outside OpenMP parallel
//regions nest: each one's name is appended to those of the timers already
//open around it, giving a path like
//    resolve/s=12/resolvor_modeled/del_and_gaussian/gaussian
//On entry it prints "DIAGBEGIN path=.. t=..", on exit one line
//    DIAG path=.. wall=.. cpu=.. threads=.. util=.. rss_gb=.. peak_gb=.. t=..
//to stderr. cpu is getrusage(RUSAGE_SELF) (all threads in this process),
//threads defaults to omp_get_max_threads(), and util = cpu/(wall*threads),
//so a fully serial step shows util ~ 1/threads. rss_gb is resident memory
//at exit, peak_gb the process peak so far, t seconds since the first timer.
//scripts/diag_summary.py aggregates these lines into a tree.
#pragma once
#include <chrono>
#include <sys/resource.h>
#include <unistd.h>
#include <cstdio>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <omp.h>

struct Diag {
	static bool enabled(){
		static bool e = std::getenv("DIAG_TIMING") != nullptr;
		return e;
	}
	static double cpu_seconds(){
		struct rusage ru;
		getrusage(RUSAGE_SELF, &ru);
		return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec*1e-6
		     + ru.ru_stime.tv_sec + ru.ru_stime.tv_usec*1e-6;
	}
	//per-calling-thread cpu, unlike cpu_seconds() above which sums the
	//whole process -- lets a parallel region attribute cpu to exactly the
	//thread that did the work, immune to other threads' idle/spin time
	static double thread_cpu_seconds(){
		struct rusage ru;
		getrusage(RUSAGE_THREAD, &ru);
		return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec*1e-6
		     + ru.ru_stime.tv_sec + ru.ru_stime.tv_usec*1e-6;
	}
	static double peak_rss_gb(){
		struct rusage ru;
		getrusage(RUSAGE_SELF, &ru);
		return ru.ru_maxrss / 1048576.0; //ru_maxrss is in KB on Linux
	}
	static double rss_gb(){
		long pages_total = 0, pages_resident = 0;
		FILE *f = std::fopen("/proc/self/statm", "r");
		if(!f) return 0;
		if(std::fscanf(f, "%ld %ld", &pages_total, &pages_resident) != 2) pages_resident = 0;
		std::fclose(f);
		return pages_resident * (double)sysconf(_SC_PAGESIZE) / (1024.0*1024.0*1024.0);
	}
	static double elapsed(){
		static auto start = std::chrono::steady_clock::now();
		return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
	}
	//names of the timers currently open outside parallel regions
	static std::vector<std::string>& stack(){
		static std::vector<std::string> s;
		return s;
	}
	static std::string path(){
		std::string p;
		for(auto const &n : stack()){
			if(!p.empty()) p += '/';
			p += n;
		}
		return p;
	}
};

struct DiagTimer {
	std::string path;
	int threads;
	bool on;
	bool nested;
	std::chrono::steady_clock::time_point t0;
	double c0;
	DiagTimer(std::string n, int nthreads=0) : threads(nthreads), on(Diag::enabled()), nested(false){
		if(!on) return;
		Diag::elapsed();
		if(threads<=0) threads = omp_get_max_threads();
		if(!omp_in_parallel()){
			Diag::stack().push_back(std::move(n));
			nested = true;
			path = Diag::path();
		}
		else
			path = n;
		std::cerr << "DIAGBEGIN path=" << path << " t=" << Diag::elapsed() << std::endl;
		t0 = std::chrono::steady_clock::now();
		c0 = Diag::cpu_seconds();
	}
	DiagTimer(DiagTimer const&) = delete;
	DiagTimer& operator=(DiagTimer const&) = delete;
	~DiagTimer(){
		if(!on) return;
		double wall = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
		double cpu = Diag::cpu_seconds()-c0;
		double util = wall>0 ? cpu/(wall*threads) : 0;
		std::cerr << "DIAG path=" << path << " wall=" << wall << " cpu=" << cpu
		          << " threads=" << threads << " util=" << util
		          << " rss_gb=" << Diag::rss_gb() << " peak_gb=" << Diag::peak_rss_gb()
		          << " t=" << Diag::elapsed() << std::endl;
		if(nested) Diag::stack().pop_back();
	}
};
