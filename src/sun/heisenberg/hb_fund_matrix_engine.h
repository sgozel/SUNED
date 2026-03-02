// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_HB_FUND_MATRIX_ENGINE_H
#define SUN_HB_FUND_MATRIX_ENGINE_H

#include "hb_fund_engine.h"

#include <vector>


namespace sun {

typedef int32_t typePk; // should use int64_t when dimension_-1 > std::numeric_limits<int32_t>::max()

class HBFundMatrixEngine : public HBFundEngine
{
public:
	HBFundMatrixEngine(nlohmann::json const& inputParam);

	void build_matrix_lookups();
	void multiply(const std::vector<double> &, std::vector<double> &, const double &, const std::string &) const override;
	
	// Multiply should be:
	//     multiply(w, u, a, method) : u <---- H*w - a*u
    
	template <class coeff_t>
    void multiply_v1_openmp(const std::vector<coeff_t>&, std::vector<coeff_t>&, const double) const;

private:
	void free_basis();
	void dump_matrix(const unsigned int k) const;
	bool load_matrix(const unsigned int k);

private:
	std::vector<std::vector<typePk>> P_;
	bool dump_matrices_;
	std::string matrix_dump_path_;
};

} // namespace sun

#endif
