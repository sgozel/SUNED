// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <chrono>
#include <stdexcept>
#include <cstdint>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
#include "common/input_parser.h"
#include "common/time.h"
#include "common/datatypes.h"
#include "sun/irrep/irrep.h"
#include "sun/utils/utils.h"

#ifdef SG_USE_VSYT
#include "sun/syt/vsyt.h"
#include "sun/syt_usage/vsyt_usage.h"
#else
#include "sun/syt/bsyt.h"
#include "sun/syt_usage/bsyt_usage.h"
#endif


int main(int argc, char* argv[])
{
	PRINT_SUNED_VERSION
	
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // READ INPUT PARAMETERS
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
	json inputParam = parseInputArguments(argc, argv);
	
    if (!inputParam.contains("alpha")) {
		throw std::runtime_error("Missing input parameter alpha. Aborting.");
	}

	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	// RUN CALCULATION
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	sun::Irrep alpha(inputParam["alpha"].get<std::vector<unsigned int>>());
    
    UINT64 falpha = sun::multiplicity(alpha);
    
    std::cout << "alpha = " << std::endl;
    alpha.print();
    std::cout << std::endl;
    std::cout << "falpha = " << falpha << std::endl;
    
	auto t0 = std::chrono::high_resolution_clock::now();
	#ifdef SG_USE_VSYT
	std::vector<sun::Int8vSYT> Y = sun::get_SYT<int8_t>(alpha);
	#else
    std::vector<sun::tbSYT> Y = sun::get_SYT(alpha);
	#endif
	time(t0, std::string("get_SYT"));
	
	/*
    for (unsigned int i = 0; i < Y.size(); ++i) {
		Y[i].print(alpha.n());
	}
	*/
    
	return 0;
}
