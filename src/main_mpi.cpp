// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
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
	
	engine.initEngine();
	engine.build_matrix_lookups();
	
	engine.eig("multiply_mpi_matrix_v1");
	
	MPI_Finalize();
	
	return 0;
}
