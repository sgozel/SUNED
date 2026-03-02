// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
#include "./common/datatypes.h"
#include "./sun/irrep/irrep.h"
#include "./sun/utils/utils.h"
#include "./sun/syt/bsyt.h"
#include "./sun/syt_usage/bsyt_usage.h"


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

	if (!inputParam.contains("alpha")) {
		throw std::runtime_error("Missing alpha in input .json file.");
	}
	sun::Irrep alpha(inputParam["alpha"].get<std::vector<unsigned int>>());
    
    UINT64 falpha = sun::multiplicity(alpha);
    
    std::cout << "alpha = " << std::endl;
    alpha.print();
    std::cout << std::endl;
    std::cout << "falpha = " << falpha << std::endl;
    
    std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
    
    for (unsigned int i=0; i<Y.size(); ++i) {
		Y[i].print(alpha.n());
	}
    
	return 0;
}
