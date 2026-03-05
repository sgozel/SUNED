// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_HB_FUND_MATRIX_ENGINE_H
#define SUN_HB_FUND_MATRIX_ENGINE_H

#include "hb_fund_engine.h"

#include <vector>

#include "../../common/numa.h"


namespace sun {

typedef int32_t typePk; // should use int64_t when dimension_-1 > std::numeric_limits<int32_t>::max()

class HBFundMatrixEngine : public HBFundEngine
{
public:
	HBFundMatrixEngine(nlohmann::json const& inputParam);

	void initEngine() override;
	void build_matrix_lookups();
	void multiply(const sg_vec<double> &, sg_vec<double> &, const double &, const std::string &) const override;
	
	// Multiply should be:
	//     multiply(w, u, a, method) : u <---- H*w - a*u
    
	template <class coeff_t>
    void multiply_v1_openmp(const sg_vec<coeff_t>&, sg_vec<coeff_t>&, const double) const;
    
    template <class coeff_t>
    void multiply_v1_openmp_numa(const sg_vec<coeff_t>&, sg_vec<coeff_t>&, const double) const;

private:
	void free_basis();
	void dump_matrix(const unsigned int k) const;
	bool load_matrix(const unsigned int k);

private:
#ifdef SG_USE_NUMA
	std::vector<sg_vec<typePk>> P_;
#else
	std::vector<std::vector<typePk>> P_;
#endif
	bool dump_matrices_;
	std::string matrix_dump_path_;
    
    mutable sg_vec<double> work_;
};

} // namespace sun

#endif
