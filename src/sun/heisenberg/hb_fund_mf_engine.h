// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_HB_FUND_MATRIX_FREE_ENGINE_H
#define SUN_HB_FUND_MATRIX_FREE_ENGINE_H

#include "hb_fund_engine.h"


namespace sun {

class HBFundMatrixFreeEngine : public HBFundEngine
{
public:
	HBFundMatrixFreeEngine(nlohmann::json const& inputParam);

	void multiply(const std::vector<double> &, std::vector<double> &, const double &, const std::string &) const override;
	
	// Multiply should be:
	//     multiply(w, u, a, method) : u <---- H*w - a*u

	template <class coeff_t>
    void multiply_v1_openmp(const std::vector<coeff_t>&, std::vector<coeff_t>&, const double) const;
    
};

} // namespace sun

#endif
