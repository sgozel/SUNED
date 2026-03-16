// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef ED_SOLVER_H
#define ED_SOLVER_H

#include <string>
#include <utility>

#include <nlohmann/json.hpp>

#include "../common/datatypes.h"
#include "../common/numa.h"
#include "../tmatrix/tmatrix.h"
#include "../lanczos/lanczosparams.h"


class EDSolver
{
public:
	EDSolver(nlohmann::json const& inputParam);

	virtual void init() = 0;
	virtual void multiply(const sg_vec<double> &, sg_vec<double> &, const double &, const std::string &) const = 0;
	// Multiply should be:
	//     multiply(w, u, a, method) : u <---- H*w - a*u
	
	double eigenvalue(const std::string &);
	
	std::pair<double, sg_vec<double>> eigenpair(const std::string &);
	
	double check_eigvec(const std::pair<double, sg_vec<double>>&, const std::string & method) const;
	
	std::pair<double, sg_vec<double>> load_eigpair(const unsigned int index) const;
	
protected:
	UINT64 dimension_;
	lanczos::LanczosParams lanczosparams_;
	
	std::string eigvec_folder_;
};

#endif
