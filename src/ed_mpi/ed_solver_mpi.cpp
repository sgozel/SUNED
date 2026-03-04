// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "ed_solver_mpi.h"

#include <iostream>
#include <iomanip>
#include <numeric>
#include <chrono>
#include <mpi.h>

#include "../lanczos_mpi/lanczos_mpi.h"


EDSolverMPI::EDSolverMPI(nlohmann::json const& inputParam)
{	
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size_);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
	
    num_threads_ = inputParam.value("num_threads", 1);
    lanczosparams_ = lanczosmpi::LanczosParams(inputParam);
    
    mpi_dimensions_.resize(mpi_world_size_);
    mpi_start_index_.resize(mpi_world_size_);
	mpi_end_index_.resize(mpi_world_size_);
}


// Returns the MPI rank to which an element of the Hilbert space belongs to
int EDSolverMPI::mpi_rank_from_index(const UINT64 i) const {
	return i/mpi_bare_dimension_;
}


// Returns the local index associated to a global index
UINT64 EDSolverMPI::mpi_local_index_from_global_index(const UINT64 i) const
{
	int rank = mpi_rank_from_index(i);
	UINT64 local_i = i - rank*mpi_bare_dimension_;
	return local_i;
}


void EDSolverMPI::mpi_get_local_dimension()
{	
	mpi_bare_dimension_ = (dimension_ + mpi_world_size_ - 1)/mpi_world_size_;
	mpi_dimension_ = (mpi_rank_==mpi_world_size_-1 ? dimension_-(mpi_world_size_-1)*mpi_bare_dimension_ : mpi_bare_dimension_);
	
	for (int rank=0; rank<mpi_world_size_; ++rank) {
		mpi_start_index_[rank] = rank*mpi_bare_dimension_;
		mpi_end_index_[rank] = (rank==mpi_world_size_-1 ? dimension_ : (rank+1)*mpi_bare_dimension_);
		mpi_dimensions_[rank] = (rank==mpi_world_size_-1 ? (dimension_-(mpi_world_size_-1)*mpi_bare_dimension_) : mpi_bare_dimension_);
	}
}


void EDSolverMPI::print_mpi_details() const
{
	if (mpi_rank_ == 0) {
		std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << "dimension = " << dimension_ << std::endl;
		std::cout << "mpi_world_size_ = " << mpi_world_size_ << std::endl;
		std::cout << "mpi_bare_dimension_ = " << mpi_bare_dimension_ << std::endl;
		std::cout << "-----------------------------------------" << std::endl;
		std::cout << "Dimensions of each rank: " << std::endl;
		for (int rank=0; rank<mpi_world_size_; ++rank) {
			std::cout << "mpi_dimension[" << rank << "] = " << mpi_dimensions_[rank] << std::endl;
		}
		std::cout << "-----------------------------------------" << std::endl;
		std::cout << "Indices per rank (start included; end excluded)" << std::endl;
		for (int rank=0; rank<mpi_world_size_; ++rank) {
			std::cout << "rank " << rank << ": " << mpi_start_index_[rank] << " ---> " << mpi_end_index_[rank] << std::endl;
		}
		std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
	}
}


double EDSolverMPI::eig(const std::string & method)
{
	if (mpi_rank_ == 0) {
		std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << "Start eig() ..." << std::endl;
		std::cout << "dimension = " << dimension_ << std::endl;
	}

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

    auto converged = [&lanczosparams=lanczosparams_](lanczosmpi::Tmatrix & tmat) -> bool {
        return lanczosmpi::convergence(tmat, lanczosparams);
    };
    
    // multiplication for 3 Lanczos vectors
    auto mult = [&method,this](const std::vector<double> &v, std::vector<double> &w) {
        multiply(v, w, 0.0, method); // w <--- H*v
    };
    
#ifdef SG_LANCZOS_EIGVEC
	
	if (mpi_rank_==0) {
		std::cout << "Multiply with 3 Lanczos vectors - eigenvalue + eigenvector - start" << std::endl;
	}
	
	std::vector<double> GS(0);
	lanczosmpi::Tmatrix tmat = lanczosmpi::lanczos_eigvec<double>(mult, converged, dimension_, GS, lanczosparams_);
	
	if (mpi_rank_==0) {
		std::cout << "Multiply with 3 Lanczos vectors - eigenvalue + eigenvector - done" << std::endl;
	}
	
#else
	
	if (mpi_rank_==0) {
		std::cout << "Multiply with 3 Lanczos vectors - eigenvalue only - start" << std::endl;
	}
	
	lanczosmpi::Tmatrix tmat = lanczosmpi::lanczos<double>(mult, converged, mpi_dimension_, lanczosparams_);
	
	if (mpi_rank_==0) {
		std::cout << "Multiply with 3 Lanczos vectors - eigenvalue only - done" << std::endl;
	}
	
#endif
	
	if (tmat.size() == 0) {
        std::cerr << "ERROR : EDSolverMPI : output Tmatrix is zero dimensional." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
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
	if (mpi_rank_ == 0) {
		std::cout << "----------------------------" << std::endl;
		std::cout << std::fixed;
		std::cout << std::setprecision(2);
		std::cout << "eig time = " << std::setw(9) << std::right << t_total << " ms" << std::endl;
	}
	
    return eigvals[0];
}


template <typename type_mult>
void EDSolverMPI::check_eigvec(lanczosmpi::Tmatrix & tmat, const std::vector<double>& GS, type_mult mult) const
{
	std::vector<double> HGS(mpi_dimension_, 0.0);
	
	mult(GS, HGS); // Hv <--- H*GS - 0.0*Hv
	
	const double energy = tmat.eigenvalues()[0];
	
	double distance_loc = std::transform_reduce(
								HGS.begin(), HGS.end(),
								GS.begin(),
								0.0,
								std::plus<>{},
								[energy](double h_el, double gs_el) -> double {
									double d = h_el - energy * gs_el;
									return d * d;});
	double distance = 0.0;
	MPI_Allreduce(&distance_loc, &distance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	distance = std::sqrt(distance);
	
	if (mpi_rank_ == 0) {
		std::ios_base::fmtflags coutflags(std::cout.flags());
		std::cout << "rank " << mpi_rank_ 
				  << ": norm(H*GS - E*GS) = " 
				  << std::setprecision(4) << std::scientific 
				  << distance << std::endl;
		std::cout.flags(coutflags);
	}
}
