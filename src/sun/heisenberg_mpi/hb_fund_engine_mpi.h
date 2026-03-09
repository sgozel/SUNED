// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_HB_FUND_ENGINE_MPI_H
#define SUN_HB_FUND_ENGINE_MPI_H

#include "../heisenberg/hb_engine.h"

#include <vector>
#include <string>

#ifdef SG_USE_BASIC_SYT
#include "../syt/vsyt.h"
#else
#include "../syt/bsyt.h"
#endif
#include "../../common/numa.h"


namespace sun {

#ifdef SG_USE_BASIC_SYT
typedef int8_t SYTel;
typedef Int8vSYT SYT;
#else
typedef tbSYT SYT;
#endif

class HBFundEngineMPI : public HBEngine
{
public:
	HBFundEngineMPI(nlohmann::json const& inputParam);

	void init() override;
	void multiply(const sg_vec<double> &, sg_vec<double> &, const double &, const std::string &) const = 0;

protected:
	std::vector<SYT> Y_;
};

} // namespace sun

#endif
