// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef ED_SOLVER_H
#define ED_SOLVER_H

#include <string>

#include <nlohmann/json.hpp>

#include "../common/datatypes.h"
#include "../tmatrix/tmatrix.h"
#include "../lanczos/lanczosparams.h"


class EDSolver
{
public:
	EDSolver(nlohmann::json const& inputParam);

	virtual void initEngine() = 0;
	virtual void multiply(const std::vector<double> &, std::vector<double> &, const double &, const std::string &) const = 0;
	// Multiply should be:
	//     multiply(w, u, a, method) : u <---- H*w - a*u
	
	double eig(const std::string &);

private:
	template <typename type_mult>
	void check_eigvec(lanczos::Tmatrix & tmat, const std::vector<double>& GS, type_mult mult3vecs) const;

protected:
	UINT64 dimension_;
	lanczos::LanczosParams lanczosparams_;
	unsigned int num_threads_;
	
	std::string eigvec_folder_;
};

#endif
