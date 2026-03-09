// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "ed_solver.h"

#include <iostream>
#include <iomanip>
#include <numeric>
#include <chrono>
#include <stdexcept>

#include "../common/time.h"
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


double EDSolver::eigenvalue(const std::string & method)
{
	std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
	std::cout << "Start eigenvalue() ..." << std::endl;
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
    auto mult = [&method,this](const sg_vec<double> &v, sg_vec<double> &w, double a) {
        multiply(v, w, a, method); // w <--- H*v - a*w
    };
	// Lanczos 2 vectors eigenvalue only
	std::cout << "Multiply with 2 Lanczos vectors - eigenvalue only - start" << std::endl;
	lanczos::Tmatrix tmat = lanczos::lanczos_two_vectors<double>(mult, converged, dimension_, lanczosparams_);
	std::cout << "Multiply with 2 Lanczos vectors - eigenvalue only - end" << std::endl;
#else
    auto mult = [&method,this](const sg_vec<double> &v, sg_vec<double> &w) {
        multiply(v, w, 0.0, method); // w <--- H*v
    };
	// Lanczos 3 vectors eigenvalue only
	std::cout << "Multiply with 3 Lanczos vectors - eigenvalue only - start" << std::endl;
	lanczos::Tmatrix tmat = lanczos::lanczos<double>(mult, converged, dimension_, lanczosparams_);
	std::cout << "Multiply with 3 Lanczos vectors - eigenvalue only - end" << std::endl;
#endif
    
	if (tmat.size() == 0) {
        throw std::runtime_error("ERROR : EDSolver : output Tmatrix is zero dimensional.");
	}
	
	std::vector<double> eigvals = tmat.eigenvalues();

	time(t0, "eigenvalue()");

    return eigvals[0];
}


std::pair<double, sg_vec<double>> EDSolver::eigenpair(const std::string & method)
{
	std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
	std::cout << "Start eigenpair() ..." << std::endl;
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
    
    sg_vec<double> GS(0);
    
#ifdef SG_LANCZOS_TWO_VECTORS
    auto mult = [&method,this](const sg_vec<double> &v, sg_vec<double> &w, double a) {
        multiply(v, w, a, method); // w <--- H*v - a*w
    };
	// Lanczos 2 vectors with eigenvector
	std::cout << "Multiply with 2 Lanczos vectors - eigenvalue + eigenvector - start" << std::endl;
	lanczos::Tmatrix tmat = lanczos::lanczos_two_vectors_eigvec(mult, converged, dimension_, GS, lanczosparams_);
	std::cout << "Multiply with 2 Lanczos vectors - eigenvalue + eigenvector - end" << std::endl;
#else
    auto mult = [&method,this](const sg_vec<double> &v, sg_vec<double> &w) {
        multiply(v, w, 0.0, method); // w <--- H*v
    };
	// Lanczos 3 vectors with eigenvector
	std::cout << "Multiply with 3 Lanczos vectors - eigenvalue + eigenvector - start" << std::endl;
	lanczos::Tmatrix tmat = lanczos::lanczos_eigvec<double>(mult, converged, dimension_, GS, lanczosparams_);
	std::cout << "Multiply with 3 Lanczos vectors - eigenvalue + eigenvector - end" << std::endl;
#endif
    
	if (tmat.size() == 0) {
        throw std::runtime_error("ERROR : EDSolver : output Tmatrix is zero dimensional.");
	}
	
	std::vector<double> eigvals = tmat.eigenvalues();

	time(t0, "eigenpair()");

    return {eigvals[0], GS};
}


double EDSolver::check_eigvec(const std::pair<double, sg_vec<double>>& eigpair, const std::string & method) const
{
	auto mult = [&method,this](const sg_vec<double> &v, sg_vec<double> &w) {
        multiply(v, w, 0.0, method); // w <--- H*v - a*w
    };
	
	const double energy = eigpair.first;
	const sg_vec<double>& GS = eigpair.second;

	sg_vec<double> HGS(dimension_);
	mult(GS, HGS); // HGS <--- H*GS

#ifdef SG_USE_NUMA
	double distance = 0.0;
	#pragma omp parallel for reduction(+:distance) schedule(static)
	for (UINT64 i = 0; i < dimension_; ++i) {
		double d = HGS[i] - energy*GS[i];
		distance += d * d;
	}
#else
	double distance = std::transform_reduce(HGS.begin(), HGS.end(), GS.begin(), 0.0,
											std::plus<>(),
											[energy](double hgsel, double gsel) -> double {
												const double d = hgsel - energy*gsel;
												return d*d;});
#endif
	distance = std::sqrt(distance);
	std::ios_base::fmtflags coutflags(std::cout.flags());
	std::cout << "norm(H*GS - E*GS) = " 
			  << std::setprecision(4) << std::scientific 
			  << distance << std::endl;
	std::cout.flags(coutflags);
	
	return distance;
}
