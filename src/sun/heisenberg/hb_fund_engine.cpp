// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "hb_fund_mf_engine.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <omp.h>

#include "../../common/time.h"
#include "../utils/utils.h"

#ifdef SG_USE_VSYT
#include "../syt_usage/vsyt_usage.h"
#else
#include "../syt_usage/bsyt_usage.h"
#endif

namespace sun {

HBFundEngine::HBFundEngine(nlohmann::json const& inputParam)
: HBEngine(inputParam)
{
	if (Ns_!=alpha_.n()) {
		throw std::runtime_error("For the fundamental irrep at each site, Ns must match the number of boxes in alpha.");
	}
	
	#ifndef SG_USE_VSYT
	tbSYT::check(N_, alpha_.n());
	#endif
	
	std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
	std::cout << "N = " << N_ << std::endl;
	std::cout << "Ns = " << Ns_ << std::endl;
	std::cout << "Target irrep: " << std::endl;
	alpha_.print();
	std::cout << std::endl;
	std::cout << "num_threads = " << omp_get_max_threads() << std::endl;
	lanczosparams_.print();
	lattice_.print_sites();
	lattice_.print_bonds();
	//lattice_.print_bonds_light();
	//lattice_.print_caching_info();
	//lattice_.print_parents();
	//lattice_.print_children();
}


void HBFundEngine::init()
{	
	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	
	dimension_ = multiplicity(alpha_);
	std::cout << "dimension = " << dimension_ << std::endl;

	#ifdef SG_USE_VSYT
	Y_ = get_SYT<SYTel>(alpha_);
	#else
	Y_ = get_SYT(alpha_);
	#endif
	
	double factor = 1e6;
	std::string F("MB");
	if (8*Y_.size()>=1e9) {
		factor *= 1000;
		std::string F("GB");
	}
	
	#ifdef SG_USE_VSYT
	std::cout << "SYTs Memory: " << alpha_.n()*sizeof(Y_[0][0])*((double)Y_.size()/factor) << F << std::endl;
	#else
	std::cout << "SYTs Memory: " << sizeof(Y_[0])*((double)Y_.size()/factor) << F << std::endl;
	#endif
	
	time(t0, std::string("init"));
}

} // namespace sun
