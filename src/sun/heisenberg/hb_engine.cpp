// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "hb_engine.h"

#include <iostream>
#include <map>
#ifdef SG_USE_MPI
#include <mpi.h>
#else
#include <stdexcept>
#endif


namespace sun {

HBEngine::HBEngine(nlohmann::json const& inputParam)
#ifdef SG_USE_MPI
: EDSolverMPI(inputParam)
#else
: EDSolver(inputParam)
#endif
{
	if (!inputParam.contains("N")) {
#ifdef SG_USE_MPI
		std::cerr << "Missing N in input .json file." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
#else
        throw std::runtime_error("Missing N in input .json file.");
#endif
    }
    N_ = inputParam["N"];
	
	if (!inputParam.contains("Ns")) {
#ifdef SG_USE_MPI
		std::cerr << "Missing Ns in input .json file." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
#else
        throw std::runtime_error("Missing Ns in input .json file.");
#endif
    }
    Ns_ = inputParam["Ns"];
	
	if (!inputParam.contains("alpha")) {
#ifdef SG_USE_MPI
		std::cerr << "Missing alpha in input .json file." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
#else
		throw std::runtime_error("Missing alpha in input .json file.");
#endif
	}
	alpha_ = Irrep(inputParam["alpha"].get<std::vector<unsigned int>>());
	if (alpha_.nrows()>N_) {
#ifdef SG_USE_MPI
		std::cerr << "The number of rows in alpha is greater than N." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
#else
		throw std::runtime_error("The number of rows in alpha is greater than N.");
#endif
	}
	
	if (!inputParam.contains("latticefile")) {
#ifdef SG_USE_MPI
		std::cerr << "Missing latticefile in input .json file." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
#else
		throw std::runtime_error("Missing latticefile in input .json file.");
#endif
	}
	std::string latticefile = inputParam["latticefile"];

	try {
		lattice_ = sun::Lattice(latticefile);
	} catch (const std::runtime_error& error) {
		std::cerr << "Reading lattice file failed with error:" << std::endl;
		std::cerr << error.what() << std::endl;
		std::cerr << "Interrupting code execution" << std::endl;
#ifdef SG_USE_MPI
		MPI_Abort(MPI_COMM_WORLD, 1);
#else
		std::abort();
#endif
	}
	
	if (Ns_!=lattice_.get_Ns()) {
#ifdef SG_USE_MPI
		std::cerr << "Problem: lattice files does not match Ns provided in .json file." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
#else
		throw std::runtime_error("Problem: lattice files does not match Ns provided in .json file.");
#endif
	}

	std::vector<std::string> couplingNames = lattice_.get_couplingNames();
	std::map<std::string, double> couplingsMap;

	for (unsigned int i=0; i<couplingNames.size(); ++i) {
		if (!inputParam.contains(couplingNames[i])) {
#ifdef SG_USE_MPI
			std::cerr << "Missing coupling in input .json file." << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
#else
			throw std::runtime_error("Missing coupling in input .json file.");
#endif
		}
		couplingsMap[couplingNames[i]] = inputParam[couplingNames[i]];
	}
	lattice_.assignCouplingValues(couplingsMap);
}

} // namespace sun
