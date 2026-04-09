// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "hb_fund_mf_engine.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>
#include <random>
#include <stdexcept>
#include <omp.h>
#include <map>
#include <unordered_map>

#include "../../common/time.h"
#include "../utils/utils.h"


#ifdef SG_USE_VSYT
#include "../syt_usage/vsyt_usage.h"
#else
#include "../syt_usage/bsyt_usage.h"
#endif


namespace sun {

HBFundMatrixFreeEngine::HBFundMatrixFreeEngine(nlohmann::json const& inputParam)
: HBFundEngine(inputParam)
{
}


void HBFundMatrixFreeEngine::multiply(const sg_vec<double> & w, sg_vec<double> & u, const double & a, const std::string & method) const
{
    if (method=="multiply_v1_openmp") {
        multiply_v1_openmp(w, u, a);
	} else if (method=="multiply_v1_openmp_numa") {
        multiply_v1_openmp_numa(w, u, a);
	} else {
        throw std::runtime_error("Multiply method undefined");
    }
}


template <class coeff_t>
void HBFundMatrixFreeEngine::multiply_v1_openmp(const sg_vec<coeff_t>& w, sg_vec<coeff_t>& u, const double a) const
{
	// u <--- H*w - a*u
	
	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	
	std::for_each(u.begin(), u.end(), [a](coeff_t& el) { el*=(-a);});

	const unsigned int fsd_freq = 2;

	UINT64 sorted_up_to = 0;

	for (size_t b = 0; b < lattice_.bonds.size(); ++b)
	{	
		//std::chrono::time_point<std::chrono::high_resolution_clock> tb0 = std::chrono::high_resolution_clock::now();
		const auto& bond = lattice_.bonds[b];
		const double J = bond.couplingValue;
		
		#pragma omp parallel for schedule(guided)
		for (UINT64 i=0; i<dimension_; ++i)
		{	
			std::vector<SYT> ydev;
			std::vector<double> coeffdev;
			ydev.push_back(Y_[i]);
			coeffdev.push_back(1.0);
			
			sorted_up_to = 0;
			bool fullsimplified = false;
			const unsigned int mid = (bond.ops.size()-1)/2;

			for (unsigned int j=0; j<bond.ops.size(); ++j) {
				int k = bond.ops[j].getk();
				develop_consecutive_number_inplace(alpha_, ydev, coeffdev, k);
				fullsimplified = false;
				if ((j>mid) && ((j+1)%fsd_freq==0)) {
					fullsimplify_development(ydev, coeffdev, sorted_up_to);
					fullsimplified = true;
					sorted_up_to = ydev.size();
				}
			}
			
			if (fullsimplified==false) {
				fullsimplify_development(ydev, coeffdev, sorted_up_to);
			}
			
			//========================
			// UPDATE LANCZOS VECTOR
			//========================
			const UINT64 Nyout = ydev.size();
			auto it = Y_.begin();
			for (UINT64 t=0; t<Nyout; ++t) {
				it = std::lower_bound(it, Y_.end(), ydev[t]);
				const UINT64 index = it - Y_.begin();
				u[i] += J * coeffdev[t] * w[index];
				++it;
			}
		}
		//time(tb0, std::string("Bond ")+std::to_string(b)+"/"+std::to_string(lattice_.get_nbonds())+": ");
	}

	time(t0, "multiply");
}



template <class coeff_t>
void HBFundMatrixFreeEngine::multiply_v1_openmp_numa(const sg_vec<coeff_t>& w, sg_vec<coeff_t>& u, const double a) const
{
	// u <--- H*w - a*u
	
	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel for schedule(static)
	for (UINT64 i = 0; i < dimension_; ++i) {
		u[i] *= -a;
	}

	const unsigned int fsd_freq = 2;

	UINT64 sorted_up_to = 0;

	for (size_t b = 0; b < lattice_.bonds.size(); ++b)
	{	
		//std::chrono::time_point<std::chrono::high_resolution_clock> tb0 = std::chrono::high_resolution_clock::now();
		const auto& bond = lattice_.bonds[b];
		const double J = bond.couplingValue;
		
		#pragma omp parallel for schedule(static)
		for (UINT64 i=0; i<dimension_; ++i)
		{	
			std::vector<SYT> ydev;
			std::vector<double> coeffdev;
			ydev.push_back(Y_[i]);
			coeffdev.push_back(1.0);
			
			sorted_up_to = 0;
			bool fullsimplified = false;
			const unsigned int mid = (bond.ops.size()-1)/2;

			for (unsigned int j=0; j<bond.ops.size(); ++j) {
				int k = bond.ops[j].getk();
				develop_consecutive_number_inplace(alpha_, ydev, coeffdev, k);
				fullsimplified = false;
				if ((j>mid) && ((j+1)%fsd_freq==0)) {
					fullsimplify_development(ydev, coeffdev, sorted_up_to);
					fullsimplified = true;
					sorted_up_to = ydev.size();
				}
			}
			
			if (fullsimplified==false) {
				fullsimplify_development(ydev, coeffdev, sorted_up_to);
			}
			
			//========================
			// UPDATE LANCZOS VECTOR
			//========================
			const UINT64 Nyout = ydev.size();
			auto it = Y_.begin();
			for (UINT64 t=0; t<Nyout; ++t) {
				it = std::lower_bound(it, Y_.end(), ydev[t]);
				const UINT64 index = it - Y_.begin();
				u[i] += J * coeffdev[t] * w[index];
				++it;
			}
		}
		//time(tb0, std::string("Bond ")+std::to_string(b)+"/"+std::to_string(lattice_.get_nbonds())+": ");
	}
	
	time(t0, "multiply");
}

} // namespace sun
