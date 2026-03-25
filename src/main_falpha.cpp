// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
#include "common/input_parser.h"
#include "common/datatypes.h"
#include "sun/irrep/irrep.h"
#include "sun/utils/utils.h"


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
    
    double factor = 1e6;
    std::string units = "MB";
    if (falpha >= 1.125e8) {
		factor *= 1000;
		units = "GB";
	}
    
    unsigned int extent = 1;
	unsigned int Nmax = 2;
	while (alpha.nrows()>Nmax) {
		extent += 1;
		Nmax = 1<<extent;
	}
    
    double lanczos_space = sizeof(double)*falpha/factor;
    double vsyt_space = alpha.n()*sizeof(int8_t)*falpha/factor;
    double bsyt_space = sizeof(UINT64)*falpha/factor;
    unsigned int bsyt_nbits = 64;
    if ((alpha.n()*extent>bsyt_nbits) && (alpha.n()*extent<=128)) {
		bsyt_nbits = 128;
		bsyt_space *= 2;
	} else if (alpha.n()*extent>128) {
		throw std::runtime_error("Huge system. Need more than 128 bits on advanced storing strategy.");
	}
    
	std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
	std::cout << "MEMORY ANALYSIS" << std::endl;
	std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
	std::cout << "1 real double-precision Lanczos vector: " << lanczos_space << units << std::endl;
	std::cout << "SYTs naive storage, int8_t: " << vsyt_space << units << std::endl;
	std::cout << "SYTs advanced storage, extent=" << extent << " : " << bsyt_space << units << std::endl;
	std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
	std::cout << "ADVANCED STORAGE STRATEGY" << std::endl;
	std::cout << "SYTs + 2 Lanczos vectors: " << 2*lanczos_space + bsyt_space << units << std::endl;
	std::cout << "SYTs + 3 Lanczos vectors: " << 3*lanczos_space + bsyt_space << units << std::endl;
	std::cout << "SYTs + 4 Lanczos vectors: " << 4*lanczos_space + bsyt_space << units << std::endl;
	std::cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::" << std::endl;
	
	return 0;
}
