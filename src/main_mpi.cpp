// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <utility> // std::pair
#include <mpi.h>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
#include "common/input_parser.h"
#include "common/numa.h"
#include "sun/heisenberg_mpi/hb_fund_matrix_engine_mpi.h"


int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // READ INPUT PARAMETERS
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
	json inputParam = parseInputArguments(argc, argv);
	
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	// RUN SIMULATION
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	sun::HBFundMatrixEngineMPI engine(inputParam);
	
	// Initialization
	engine.init();
	engine.build_matrix_lookups();

	const std::string mvm_method("multiply_mpi_matrix_v1");
	
	// Compute eigenvalues
	double energy = engine.eigenvalue(mvm_method);
	if (mpi_rank == 0) {
		std::cout << "In main: from eigenvalue(" << mvm_method << "), energy = " << energy << std::endl;
	}
	
	
	// Compute eigenpair
	std::pair<double, sg_vec<double>> eigpair = engine.eigenpair(mvm_method);
    if (mpi_rank == 0) {
		std::cout << "In main: from eigenpair(), energy = " << eigpair.first << std::endl;
	}
    engine.check_eigvec(eigpair, mvm_method);
    
	// Compute correlations
	const unsigned int refsite = inputParam.value("correlation_refsite", 0);
    engine.correlations(eigpair.second, refsite);
	
	if (mpi_rank == 0) {
		std::cout << "Leaving main" << std::endl;
	}
	
	MPI_Finalize();
	
	return 0;
}
