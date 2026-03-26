// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "tmatrix.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>


extern "C" {
void dstev_(const char* JOBZ, // 'N': eigenvalues only; 'V': eigenvalues and eigenvectors.
			const int* N,	  // matrix size
			double* D,	      // size N: on entry: diagonal values of tridiagonal matrix; on exit, if INFO = 0, eigenvalues in ascending order
			double* E,	      // size N-1: on entry, subdiagonal elements of the tridiagonal matrix
			double* Z,	      // If JOBZ = 'V': on output: orthonormal eigenvectors of the matrix, in flat format
			const int* LDZ,   // Leading dimension of the array Z. LDZ > = 1, and if JOBZ = 'V', LDZ > = max(1, N)
			double* WORK,	  // WORK dimension max(1, 2*N-2)
			int* INFO);	      // output status: = 0: successful exit; < 0: if INFO = -i, i-th argument had illegal value; > 0: if INFO  = i, algorithm 
							  // failed to converge; i off-diagonal elements of OFFD did not converge to zero
}

#ifdef SG_USE_MPI
namespace lanczosmpi {
#else
namespace lanczos {
#endif

Tmatrix::Tmatrix()
: n_(0)
{
}

Tmatrix::Tmatrix(const std::vector<double>& alpha, const std::vector<double> & beta)
: alpha_(alpha), beta_(beta), n_(alpha.size())
{
    if (n_!=beta_.size()) {
		throw std::runtime_error("Tmatrix: mismatched sizes");
    }
}

void Tmatrix::push_back(const double& a, const double & b)
{
    alpha_.push_back(a);
    beta_.push_back(b);
    n_ += 1;
}

void Tmatrix::pop_back()
{
    alpha_.pop_back();
    beta_.pop_back();
    n_ -= 1;
}

std::vector<double> Tmatrix::eigenvalues()
{
	std::vector<double> eigvals(alpha_);
	
	if (eigvals_.count(n_)) {
		eigvals = eigvals_[n_];
	} else {
		if (n_==0) {
			throw std::runtime_error("tmat.eigenvalues() : tmat is zero-dimensional.");
		} else if (n_>1) {
			
			const char JOBZ = 'N'; // eigenvalues
			int N = n_;
			int LDZ = 1;
			int INFO = 0;
			
			std::vector<double> subdiag(beta_);
			subdiag.pop_back();
			
			double* Z = nullptr;
			std::vector<double> WORK(std::max(static_cast<unsigned int>(1), 2*n_ - 2));
			
			dstev_(&JOBZ,
				   &N, 
				   eigvals.data(),	
				   subdiag.data(), 
				   Z,   
				   &LDZ, 
				   WORK.data(),	
				   &INFO);
			
			if (INFO != 0) {
				throw std::runtime_error("dstev_ failed with INFO = " + std::to_string(INFO));
			}
		}
		eigvals_[n_] = eigvals;
	}
	
	return eigvals;
}

std::vector<std::vector<double>> Tmatrix::eigenvectors()
{
    return eig().second;
}

std::pair<std::vector<double>, std::vector<std::vector<double>>> Tmatrix::eig()
{
	std::vector<double> eigvals(alpha_);
	std::vector<std::vector<double>> eigvecs(n_, std::vector<double>(n_));
	
	if (n_==0) {
		throw std::runtime_error("tmat.eigenvectors() : tmat is zero-dimensional.");
	} else if (n_==1) {
		eigvecs[0][0] = 1.0;
	} else {
		
		const char JOBZ = 'V'; // eigenvectors
		int N = n_;
		int LDZ = N;
		int INFO = 0;
		
		std::vector<double> subdiag(beta_);
		subdiag.pop_back();
		
		std::vector<double> Z(n_ * n_); // eigenvectors stored column-wise, flat format
		std::vector<double> WORK(std::max(static_cast<unsigned int>(1), 2*n_ - 2));
		
		dstev_(&JOBZ, 
			   &N, 
			   eigvals.data(), 
			   subdiag.data(), 
			   Z.data(), 
			   &LDZ, 
			   WORK.data(),
			   &INFO);
		
		if (INFO != 0) {
			throw std::runtime_error("dstev_ failed with INFO = " + std::to_string(INFO));
		}
		
		// convert to vector of eigenvectors
		for (unsigned int i = 0; i < n_; i++) {
			for (unsigned int j = 0; j < n_; j++) {
				eigvecs[i][j] = Z[i*n_+j];
			}
		}
	}
	
	eigvals_[n_] = eigvals;
    
    return {eigvals, eigvecs};
}

void Tmatrix::log(const std::string& folder) const
{
	std::string filename = folder + "tmat.log";
	
	std::ofstream out(filename);
    if (!out) {
		throw std::runtime_error("Cannot open file: " + filename);
    }
    
    out << std::fixed << std::setprecision(16);
    
    out << "[n]" << std::endl << n_ << std::endl;
    out << "[alpha]" << std::endl;
    for (unsigned int i = 0; i < n_; ++i) {
        out << alpha_[i] << std::endl;
    }
    out << "[beta]" << std::endl;
    for (unsigned int i = 0; i < n_; ++i) {
        out << beta_[i] << std::endl;
    }
}

void Tmatrix::log_eigvals(unsigned int k, const std::string& folder) const
{
	k = std::min(k, n_);
	
	for (unsigned int t = 0; t < k; ++t)
	{
		std::string eigval_filename = folder + "eigval_" + std::to_string(t) + ".log";
		std::ofstream out(eigval_filename);
		if (!out) {
			throw std::runtime_error("Cannot open file: " + eigval_filename);
		}
		out << "Eigenvalue " << t << " of Tmatrix" << std::endl;
		out << std::fixed << std::setprecision(16);
		
		for (const auto& [key, value] : eigvals_) {
			if (t < value.size()) {
				out << key << ": " << value[t] << std::endl;
			}
		}
	}
}

}
