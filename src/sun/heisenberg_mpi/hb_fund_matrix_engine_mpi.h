// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_HB_FUND_MATRIX_ENGINE_MPI_H
#define SUN_HB_FUND_MATRIX_ENGINE_MPI_H

#include "hb_fund_engine_mpi.h"

#include <vector>
#include <string>
#include <utility>

#include "../../common/numa.h"


namespace sun {

typedef int64_t typePk; // can potentially use int32_t

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
    
    template <class coeff_t>
    void multiply_mpi_matrix_v1_numa(const sg_vec<coeff_t>&, sg_vec<coeff_t>&, const double) const;

private:
	void free_basis();
	void dump_matrix(const unsigned int k) const;
	bool load_matrix(const unsigned int k);

protected:
	bool dump_matrices_;
	std::string matrix_dump_path_;
	
	std::vector<std::vector<typePk>> P_;
	
	std::vector<UINT64> mpi_nb_offdiag_;
	std::vector<std::vector<INT64>> mpi_offdiag_nodes_;
	std::vector<std::vector<INT64>> mpi_offdiag_nodes_acc_;
	std::vector<std::vector<UINT64>> mpi_local_index_base_;
	std::vector<std::vector<UINT64>> mpi_local_index_friend_;

	std::vector<std::pair<SYT, SYT>> Y_bounds_;

private:
	mutable sg_vec<double> work_;
};

} // namespace sun

#endif
