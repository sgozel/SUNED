// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "hb_fund_matrix_engine.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <chrono>
#include <stdexcept>
#include <omp.h>

#include "../../common/time.h"
#include "../utils/utils.h"
#include "../utils/young_factor.h"


#ifdef SG_USE_VSYT
#include "../syt_usage/vsyt_usage.h"
#else
#include "../syt_usage/bsyt_usage.h"
#endif


namespace sun {

HBFundMatrixEngine::HBFundMatrixEngine(nlohmann::json const& inputParam)
: HBFundEngine(inputParam)
{	
	P_.resize(alpha_.n()-1);
	
	dump_matrices_ = inputParam.value("dump_matrices", false);
	if (dump_matrices_ == true) {
		if (!inputParam.contains("matrix_dump_folder_path")) {
			throw std::runtime_error("Missing matrix_dump_folder_path in input .json file for HBFundMatrixEngine.");
		}
		matrix_dump_path_ = inputParam["matrix_dump_folder_path"];
		if (matrix_dump_path_.back()!='/') {
			matrix_dump_path_ += std::string("/");
		}
	}
}

void HBFundMatrixEngine::init()
{
	HBFundEngine::init();
	
	work_.resize(dimension_);
	#pragma omp parallel for schedule(static)
	for (UINT64 i = 0; i < dimension_; ++i) {
		work_[i] = 0.0;
	}
}


void HBFundMatrixEngine::build_matrix_lookups()
{
	std::chrono::time_point<std::chrono::high_resolution_clock> t_start = std::chrono::high_resolution_clock::now();
	
	{
		double factor = 1e6;
		std::string units = "MB";
		if (8*dimension_ >= 1e9) {
			factor *= 1000;
			units = "GB";
		}
		double memMatrixLookups = sizeof(typePk)*(alpha_.n()-1)*static_cast<double>(dimension_)/factor;
		double memLanczos = sizeof(double)*3*static_cast<double>(dimension_)/factor;
		double memLanczosCopy = sizeof(double)*static_cast<double>(dimension_)/factor;
		double memTotal = memMatrixLookups + memLanczos + memLanczosCopy;
	
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << ":::::::: TOTAL MEMORY REQUIREMENTS ::::::::" << std::endl;
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << "Matrix lookups: " << memMatrixLookups << units << std::endl;
		std::cout << "3 Lanczos vectors: " << memLanczos << units << std::endl;
		std::cout << "1 Lanczos work in multiply: " << memLanczosCopy << units << std::endl;
		std::cout << "--------------" << std::endl;
		std::cout << "Total: " << memTotal << units << std::endl;
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
	}
	
	for (unsigned int k=0; k<alpha_.n()-1; ++k)
	{	
		std::chrono::time_point<std::chrono::high_resolution_clock> tk_start = std::chrono::high_resolution_clock::now();
		
		std::string transpo_string = std::string("(") + std::to_string(k) + ", " + std::to_string(k+1) + ")";
		
		P_[k].resize(dimension_);
		
		// Attempt to load P_[k] from file if it exists
        if ((dump_matrices_ == true) && (load_matrix(k))) {
			time(tk_start, std::string("Loaded from file ")+transpo_string);
            continue;
        }
		
		//#pragma omp parallel for schedule(guided)
		#pragma omp parallel for schedule(static)
		for (UINT64 i=0; i<dimension_; ++i) {
			
			const int rowk = Y_[i].get(k);
			const int rowkk = Y_[i].get(k+1);
			
			if (rowk == rowkk) {
#ifdef SG_STORE_COLUMNS
				P_[k][i] = -1;
#else
				P_[k][i] = 0;
#endif
			} else {
				const std::pair<int, int> cy = get_column_k_k_plus_one(Y_[i], k);
				if (cy.first == cy.second) {
#ifdef SG_STORE_COLUMNS
					P_[k][i] = 0;
#else
					P_[k][i] = -1;
#endif
				} else {
					SYT yfriend = Y_[i];
					yfriend.exchange(k, k+1);
					
					auto it = std::lower_bound(Y_.begin(), Y_.end(), yfriend);
					UINT64 index = it - Y_.begin();
					
					if (index > i) {
						P_[k][i] = static_cast<typePk>(index);
						const typePk ax = cy.first - rowk - cy.second + rowkk; // axial distance from k to k+1 in SYT Y_[i]
						if (ax >= 0) {
							// This should never happen, because SYTs in Y_
							// are ordered in the descending order of the LLOS
							throw std::runtime_error("Problem: ax >= 0");
						}
						P_[k][index] = ax;
						// negative and different from -1
					}
				}
			}
		}
		
		time(tk_start, std::string("Building ")+transpo_string);
		
		if (dump_matrices_ == true) {
			dump_matrix(k);
		}
	}
	
	time(t_start, "build_matrix_lookups");
	
	free_basis();
}


void HBFundMatrixEngine::dump_matrix(const unsigned int k) const
{
	const std::string filename = matrix_dump_path_ + std::string("Pk_") + std::to_string(k) + ".bin";
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("Cannot open file : " + filename);
	}
	
	// Write dimension_
	out.write(reinterpret_cast<const char*>(&dimension_), sizeof(UINT64));
	
	// Write P_[k]
	out.write(reinterpret_cast<const char*>(P_[k].data()), dimension_ * sizeof(typePk));
}


bool HBFundMatrixEngine::load_matrix(const unsigned int k)
{
    const std::string filename = matrix_dump_path_ + std::string("Pk_") + std::to_string(k) + ".bin";
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        return false;
    }

    // Read and verify dimension_
    UINT64 stored_dimension = 0;
    in.read(reinterpret_cast<char*>(&stored_dimension), sizeof(UINT64));
    if (!in || stored_dimension != dimension_) {
        std::cerr << "Warning: dimension mismatch or read error in file: " << filename
                  << " (stored=" << stored_dimension << ", expected=" << dimension_ << ")" << std::endl;
        return false;
    }
	
    // Read P_[k]
    in.read(reinterpret_cast<char*>(P_[k].data()), dimension_ * sizeof(typePk));
    if (!in) {
        std::cerr << "Warning: failed to read matrix data from file: " << filename << std::endl;
        return false;
    }
#ifdef SG_USE_NUMA
    // Re-touch to restore NUMA layout
	{
		sg_vec<typePk> tmp(dimension_);
		#pragma omp parallel for schedule(static)
		for (UINT64 i = 0; i < dimension_; ++i) {
			tmp[i] = P_[k][i];
		}
		std::swap(P_[k], tmp);
	}
#endif
	
    return true;
}


void HBFundMatrixEngine::free_basis()
{
	std::vector<SYT>().swap(Y_);
}


void HBFundMatrixEngine::multiply(const sg_vec<double> & w, sg_vec<double> & u, const double & a, const std::string & method) const
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
void HBFundMatrixEngine::multiply_v1_openmp(const sg_vec<coeff_t>& w, sg_vec<coeff_t>& u, const double a) const
{
	// u <--- H*w - a*u
	
	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	
	std::for_each(u.begin(), u.end(), [a](coeff_t& el) {el*=(-a);});

	std::vector<coeff_t> work(dimension_); // !!! COPY OF LANCZOS VECTOR !!!

	for (size_t b = 0; b < lattice_.bonds.size(); ++b)
	{	
		//std::chrono::time_point<std::chrono::high_resolution_clock> tb0 = std::chrono::high_resolution_clock::now();
		const auto& bond = lattice_.bonds[b];
		const double J = bond.couplingValue;
		
		std::copy(w.begin(), w.end(), work.begin());
		
		for (unsigned int j=0; j<bond.ops.size(); ++j)
		{
			const unsigned int k = bond.ops[j].getk();
		
			#pragma omp parallel for schedule(guided)
			for (UINT64 i=0; i<dimension_; ++i)
			{	
				if (P_[k][i] == -1) {
					work[i] *= -1;
				} else if (P_[k][i] > 0) {
					const UINT64 index = P_[k][i];
					const double rho = 1.0/static_cast<double>(YOUNG_FACTOR * P_[k][index]);
					const double work_i = work[i];
					const double work_index = work[index];
					const double eta = std::sqrt(1.0 - rho*rho);
					work[i] = -rho * work_i + eta * work_index;
					work[index] = eta * work_i + rho * work_index;
				}
			}
		}
		
		// update Lanczos vector
		for (UINT64 i=0; i<dimension_; ++i) {
			u[i] += J*work[i];
		}
		
	}

	time(t0, "multiply");
}


template <class coeff_t>
void HBFundMatrixEngine::multiply_v1_openmp_numa(const sg_vec<coeff_t>& w, sg_vec<coeff_t>& u, const double a) const
{
	// u <--- H*w - a*u
	
	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel for schedule(static)
	for (UINT64 i = 0; i < dimension_; ++i) {
		u[i] *= -a;
	}

	for (size_t b = 0; b < lattice_.bonds.size(); ++b)
	{	
		//std::chrono::time_point<std::chrono::high_resolution_clock> tb0 = std::chrono::high_resolution_clock::now();
		const auto& bond = lattice_.bonds[b];
		const double J = bond.couplingValue;
		
		#pragma omp parallel for schedule(static)
		for (UINT64 i = 0; i < dimension_; ++i) {
			work_[i] = w[i];
		}
		
		for (unsigned int j=0; j<bond.ops.size(); ++j)
		{
			const unsigned int k = bond.ops[j].getk();
		
			#pragma omp parallel for schedule(static)
			for (UINT64 i=0; i<dimension_; ++i)
			{	
				if (P_[k][i] == -1) {
					work_[i] *= -1;
				} else if (P_[k][i] > 0) {
					const UINT64 index = P_[k][i];
					const double rho = 1.0/static_cast<double>(YOUNG_FACTOR * P_[k][index]);
					const double work_i = work_[i];
					const double work_index = work_[index];
					const double eta = std::sqrt(1.0 - rho*rho);
					work_[i] = -rho * work_i + eta * work_index;
					work_[index] = eta * work_i + rho * work_index;
				}
			}
		}
		
		// update Lanczos vector
		#pragma omp parallel for schedule(static)
		for (UINT64 i = 0; i < dimension_; ++i) {
			u[i] += J * work_[i];
		}
		
		//time(tb0, std::string("Bond ")+std::to_string(b)+"/"+std::to_string(lattice_.get_nbonds())+": ");
	}

	time(t0, "multiply");
}



} // namespace sun
