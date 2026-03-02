// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <iomanip>
#include <fstream>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
#include "sun/irrep/irrep.h"
//#include "sun/heisenberg/hb_fund_mf_engine.h"
#include "sun/heisenberg/hb_fund_matrix_engine.h"


int main(int argc, char* argv[])
{
	PRINT_SUNED_VERSION
	
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

	sun::HBFundMatrixEngine engine(inputParam);
    engine.initEngine();
    engine.build_matrix_lookups();
    engine.eig("multiply_v1_openmp");
	
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	std::cout << "leaving main" << std::endl;
	
	return 0;
}
