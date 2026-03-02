// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "hb_fund_engine_mpi.h"

#include <iostream>
#include <mpi.h>

#include "../utils/utils.h"

#ifdef SG_USE_BASIC_SYT
#include "../syt_usage/vsyt_usage.h"
#else
#include "../syt_usage/bsyt_usage.h"
#endif

namespace sun {

HBFundEngineMPI::HBFundEngineMPI(nlohmann::json const& inputParam)
: HBEngine(inputParam)
{
	if (Ns_!=alpha_.n()) {
		std::cerr << "ERROR: For the fundamental irrep at each site, Ns must match the number of boxes in alpha." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	#ifndef SG_USE_BASIC_SYT
	tbSYT::check(N_, alpha_.n());
	#endif
}


void HBFundEngineMPI::initEngine()
{	
	dimension_ = multiplicity(alpha_);
	
	if (mpi_rank_==0) {
		std::cout << "dimension = " << dimension_ << std::endl;
	}
	
	#ifdef SG_USE_BASIC_SYT
	Y_ = get_SYT<SYTel>(alpha_);
	#else
	Y_ = get_SYT(alpha_);
	#endif
	
	mpi_get_local_dimension();
	print_mpi_details();
	
	if (mpi_rank_==0) {
		double factor = 1e6;
		std::string F("MB");
		if (Y_.size()>1e9) {
			factor *= 1000;
			F = std::string("GB");
		}
		
		#ifdef SG_USE_BASIC_SYT
		std::cout << "SYTs Memory: " << alpha_.n()*sizeof(Y_[0][0])*((double)Y_.size()/factor) << F << std::endl;
		#else
		std::cout << "SYTs Memory: " << sizeof(Y_[0])*((double)Y_.size()/factor) << F << std::endl;
		#endif
	}
}

} // namespace sun
