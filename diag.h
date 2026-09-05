//diag.h -- opt-in wall/cpu timing, off by default (DIAG_TIMING env var).
//Prints "DIAG <name> wall=.. cpu=.. [threads=..]" to stderr on scope exit.
//cpu is getrusage(RUSAGE_SELF) (sums all threads in this process) delta,
//so wall*threads vs cpu gives the utilization ratio for a parallel region.
#pragma once
#include <chrono>
#include <sys/resource.h>
#include <iostream>
#include <cstdlib>
#include <string>

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
};

struct DiagTimer {
	std::string name;
	int threads;
	bool on;
	std::chrono::steady_clock::time_point t0;
	double c0;
	DiagTimer(std::string n, int nthreads=0) : name(std::move(n)), threads(nthreads), on(Diag::enabled()){
		if(on){ t0 = std::chrono::steady_clock::now(); c0 = Diag::cpu_seconds(); }
	}
	~DiagTimer(){
		if(!on) return;
		double wall = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
		double cpu = Diag::cpu_seconds()-c0;
		std::cerr << "DIAG " << name << " wall=" << wall << " cpu=" << cpu;
		if(threads>0) std::cerr << " threads=" << threads;
		std::cerr << std::endl;
	}
};
