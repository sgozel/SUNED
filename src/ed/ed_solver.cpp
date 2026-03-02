// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "ed_solver.h"

#include <iostream>
#include <iomanip>
#include <numeric>
#include <chrono>
#include <stdexcept>

#include "../lanczos/lanczos.h"


EDSolver::EDSolver(nlohmann::json const& inputParam)
{
	num_threads_ = inputParam.value("num_threads", 1);	
	lanczosparams_ = lanczos::LanczosParams(inputParam);
	eigvec_folder_ = inputParam.value("eigvec_folder_path", std::string(""));
	if ((!eigvec_folder_.empty()) && (eigvec_folder_.back()!='/')) {
		eigvec_folder_ += '/';
	}
}


double EDSolver::eig(const std::string & method)
{
	std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
	std::cout << "Start eig() ..." << std::endl;
    std::cout << "dimension = " << dimension_ << std::endl;

	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();

    if (dimension_<lanczosparams_.k) {
        lanczosparams_.k = (unsigned int) dimension_;
    }
    unsigned int nshow = 4;
    if (dimension_<nshow) {
        nshow = static_cast<unsigned int>(dimension_);
    } else {
        nshow = std::max(nshow, lanczosparams_.k);
    }

    auto converged = [&lanczosparams=lanczosparams_](lanczos::Tmatrix & tmat) -> bool {
        return lanczos::convergence(tmat, lanczosparams);
    };
    
    //:::::::::::::::::::::::::::::::::::::
    //:::::::::: DIAGONALIZATION ::::::::::
    //:::::::::::::::::::::::::::::::::::::
    
#ifdef SG_LANCZOS_TWO_VECTORS
    auto mult = [&method,this](const std::vector<double> &v, std::vector<double> &w, double a) {
        multiply(v, w, a, method); // w <--- H*v - a*w
    };
#ifdef SG_LANCZOS_EIGVEC
	// Lanczos 2 vectors with eigenvector
	std::cout << "Multiply with 2 Lanczos vectors - eigenvalue + eigenvector - start" << std::endl;
	std::vector<double> GS(0);
	lanczos::Tmatrix tmat = lanczos::lanczos_two_vectors_eigvec(mult, converged, dimension_, GS, lanczosparams_);
	bool control_eigvec = true;
	std::cout << "Multiply with 2 Lanczos vectors - eigenvalue + eigenvector - end" << std::endl;
#else
	// Lanczos 2 vectors eigenvalue only
	std::cout << "Multiply with 2 Lanczos vectors - eigenvalue only - start" << std::endl;
	lanczos::Tmatrix tmat = lanczos::lanczos_two_vectors<double>(mult, converged, dimension_, lanczosparams_);
	std::cout << "Multiply with 2 Lanczos vectors - eigenvalue only - end" << std::endl;
#endif
#else
    auto mult = [&method,this](const std::vector<double> &v, std::vector<double> &w) {
        multiply(v, w, 0.0, method); // w <--- H*v
    };
#ifdef SG_LANCZOS_EIGVEC
	// Lanczos 3 vectors with eigenvector
	std::cout << "Multiply with 3 Lanczos vectors - eigenvalue + eigenvector - start" << std::endl;
	std::vector<double> GS(0);
	lanczos::Tmatrix tmat = lanczos::lanczos_eigvec<double>(mult, converged, dimension_, GS, lanczosparams_);
	bool control_eigvec = true;
	std::cout << "Multiply with 3 Lanczos vectors - eigenvalue + eigenvector - end" << std::endl;
#else
	// Lanczos 3 vectors eigenvalue only
	std::cout << "Multiply with 3 Lanczos vectors - eigenvalue only - start" << std::endl;
	lanczos::Tmatrix tmat = lanczos::lanczos<double>(mult, converged, dimension_, lanczosparams_);
	std::cout << "Multiply with 3 Lanczos vectors - eigenvalue only - end" << std::endl;
#endif
#endif
    
	if (tmat.size() == 0) {
        throw std::runtime_error("ERROR : EDSolver : output Tmatrix is zero dimensional.");
	}
	
	std::vector<double> eigvals = tmat.eigenvalues();
	
#ifdef SG_LANCZOS_EIGVEC
	if (control_eigvec) {
		check_eigvec(tmat, GS, mult);
	}
#endif

	std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> dt_total = t1 - t0;
    double t_total = dt_total.count();
	std::cout << "----------------------------" << std::endl;
	std::cout << std::fixed;
    std::cout << std::setprecision(2);
    std::cout << "eig time = " << std::setw(9) << std::right << t_total << " ms" << std::endl;

    return eigvals[0];
}


template <typename type_mult>
void EDSolver::check_eigvec(lanczos::Tmatrix & tmat, const std::vector<double>& GS, type_mult mult) const
{
	std::vector<double> HGS(dimension_, 0.0);
	
#ifdef SG_LANCZOS_TWO_VECTORS
	mult(GS, HGS, 0.0);
#else
	mult(GS, HGS);
#endif
	
	const double energy = tmat.eigenvalues()[0];
	
	double distance = std::transform_reduce(HGS.begin(), HGS.end(), GS.begin(), 0.0,
											std::plus<>(),
											[energy](double hgsel, double gsel) -> double {
												const double d = hgsel - energy*gsel;
												return d*d;});
	distance = std::sqrt(distance);
	std::ios_base::fmtflags coutflags(std::cout.flags());
	std::cout << "norm(H*GS - E*GS) = " 
			  << std::setprecision(4) << std::scientific 
			  << distance << std::endl;
	std::cout.flags(coutflags);
}
