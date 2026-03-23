// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>
#include <mpi.h>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
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
	
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // READ INPUT PARAMETERS
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    std::ifstream inputFileStream;
    json inputParam;

    try {
        if(argc==2) {
            inputFileStream.open(argv[1]);
            inputFileStream >> inputParam;
            inputFileStream.close();
        } else {
            std::cerr << "No input file provided. Aborting." << std::endl;
            abort();
        }
    }
    catch (...) {
        std::cerr << "Caught exception at top level in [main]." << std::endl;
        abort();
    }

	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	sun::HBFundMatrixEngineMPI engine(inputParam);
	
	engine.init();
	engine.build_matrix_lookups();

	const std::string mvm_method("multiply_mpi_matrix_v1");
	
	// Compute eigenvalue
	double energy = engine.eigenvalue(mvm_method);
	if (mpi_rank == 0) {
		std::cout << "In main: from eigenvalue(" << mvm_method << "), energy = " << energy << std::endl;
	}
	
	/*
	// Compute eigenpair
	std::pair<double, sg_vec<double>> eigpair = engine.eigenpair(mvm_method);
    if (mpi_rank == 0) {
		std::cout << "In main: from eigenpair(), energy = " << eigpair.first << std::endl;
	}
    engine.check_eigvec(eigpair, mvm_method);
    
    engine.correlations(eigpair.second);
    
	*/
	
	if (mpi_rank == 0) {
		std::cout << "Leaving main" << std::endl;
	}
	
	MPI_Finalize();
	
	return 0;
}
