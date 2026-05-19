// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_HB_FUND_MATRIX_ENGINE_MPI_H
#define SUN_HB_FUND_MATRIX_ENGINE_MPI_H

#include "hb_fund_engine_mpi.h"

#include <vector>
#include <string>
#include <utility>

#include "../../common/numa.h"


namespace sun {

typedef int8_t typePk;
typedef uint32_t typeIndex;

class HBFundMatrixEngineMPI : public HBFundEngineMPI
{
public:
	HBFundMatrixEngineMPI(nlohmann::json const& inputParam);
	
	void init() override;
	void build_matrix_lookups();
	void multiply(const sg_vec<double> &, sg_vec<double> &, const double &, const std::string &) const override;
	
	// Multiply should be:
	//     multiply(w, u, a, method) : u <---- H*w - a*u
    
	template <class coeff_t>
    void multiply_mpi_matrix_v1(const sg_vec<coeff_t>&, sg_vec<coeff_t>&, const double) const;
	
	template <class T>
	std::vector<double> correlations(const T&, const unsigned int refsite = 0) const;
	
private:
	void communicate_bounds();
	void free_basis();
	void dump_matrix(const unsigned int k) const;
	bool load_matrix(const unsigned int k);
	void precise_memory_usage() const;
	
	void apply_transpositions(const std::vector<AdjacentTransposition>&) const;

protected:
	bool dump_matrices_;
	std::string matrix_dump_path_;
	
	std::vector<std::pair<SYT, SYT>> Y_bounds_;
	
	std::vector<std::vector<typePk>> P_;
	
	std::vector<UINT64> local_pairs_; // [k] count of local pairs in transposition (k, k+1)
	std::vector<UINT64> remote_pairs_; // [k] count of remote pairs in transposition (k, k+1)
	
	std::vector<UINT64> mpi_nb_offdiag_; // [k] number of off-diag elements (local + remote)
	std::vector<UINT64> mpi_nb_offdiag_local_; // [k] number of local off-diag elements
	std::vector<UINT64> mpi_nb_offdiag_remote_; // [k] number of remote off-diag elements
	
	std::vector<std::vector<INT64>> mpi_offdiag_nodes_remote_only_; // [k][rank] off-diag counts going to each rank, omitting local pairs
	std::vector<std::vector<INT64>> mpi_offdiag_nodes_remote_only_acc_; // [k][rank] accumulated off-diag counts going to each rank (displacements), omitting local pairs
	
	std::vector<std::vector<typeIndex>> mpi_local_index_base_local_;    // [k][i] local base index for local pairs
	std::vector<std::vector<typeIndex>> mpi_local_index_friend_local_;  // [k][i] local friend index for local pairs (same rank)
	
	std::vector<std::vector<typeIndex>> mpi_local_index_base_remote_;    // [k][i] local base index for remote pairs
	std::vector<std::vector<typeIndex>> mpi_local_index_friend_remote_;  // [k][i] local friend index for remote pairs (different rank)
	
	UINT64 max_offdiag_;
	UINT64 max_offdiag_remote_;
	
private:
	mutable sg_vec<double> work_;
	mutable std::vector<double> buffer_;
	mutable std::vector<double> recv_coeffs_;
};

} // namespace sun

#endif
