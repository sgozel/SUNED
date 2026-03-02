// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_HB_ENGINE_H
#define SUN_HB_ENGINE_H

#include <vector>
#include <string>

#include "../../nlohmann/json.hpp"

#ifdef SG_USE_MPI
#include "../../ed_mpi/ed_solver_mpi.h"
#else
#include "../../ed/ed_solver.h"
#endif

#include "../irrep/irrep.h"
#include "../lattice/lattice.h"


namespace sun {

#ifdef SG_USE_MPI
class HBEngine : public EDSolverMPI
#else
class HBEngine : public EDSolver
#endif
{
public:
	HBEngine(nlohmann::json const& inputParam);
	
	void initEngine() = 0;
	void multiply(const std::vector<double> &, std::vector<double> &, const double &, const std::string &) const = 0;

protected:
	unsigned int N_;
	unsigned int Ns_;
	Irrep alpha_;
	Lattice lattice_;
};

}

#endif
