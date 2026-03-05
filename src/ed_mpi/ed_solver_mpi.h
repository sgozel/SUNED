// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef ED_SOLVER_MPI_H
#define ED_SOLVER_MPI_H

#include <vector>
#include <string>
#include <utility>

#include "../nlohmann/json.hpp"

#include "../common/datatypes.h"
#include "../common/numa.h"
#include "../tmatrix/tmatrix.h"
#include "../lanczos/lanczosparams.h"


class EDSolverMPI
{
public:
	EDSolverMPI(nlohmann::json const& inputParam);

	virtual void initEngine() = 0;
	virtual void multiply(const sg_vec<double> &, sg_vec<double> &, const double &, const std::string &) const = 0;
	// Multiply should be:
	//     multiply(w, u, a, method) : u <---- H*w - a*u
	
	double eig(const std::string &);

protected:
	void mpi_get_local_dimension();
	int mpi_rank_from_index(const UINT64 i) const;
	UINT64 mpi_local_index_from_global_index(const UINT64 i) const;
	
	virtual void print_mpi_details() const;

private:
	template <typename type_mult>
	void check_eigvec(lanczosmpi::Tmatrix & tmat, const sg_vec<double>& GS, type_mult mult) const;

protected:
	UINT64 dimension_;
	lanczosmpi::LanczosParams lanczosparams_;
	unsigned int num_threads_;
	
	// MPI related attributes
	int mpi_world_size_;
	int mpi_rank_;
	UINT64 mpi_bare_dimension_;
	UINT64 mpi_dimension_;
	std::vector<UINT64> mpi_dimensions_;
	std::vector<UINT64> mpi_start_index_;
	std::vector<UINT64> mpi_end_index_;
};

#endif
