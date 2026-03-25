// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
#include "common/input_parser.h"
#include "sun/lattice/lattice.h"


int main(int argc, char* argv[])
{
	PRINT_SUNED_VERSION
	
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // READ INPUT PARAMETERS
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
	json inputParam = parseInputArguments(argc, argv);
	
    if (!inputParam.contains("latticefile")) {
		throw std::runtime_error("Missing input parameter latticefile. Aborting.");
	}
    
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	// RUN CALCULATION
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
    std::string latticefile = inputParam["latticefile"];

    sun::Lattice lattice;
    
    try {
		lattice = sun::Lattice(latticefile);
	} catch (const std::runtime_error& error) {
		std::cerr << "Reading lattice file failed with error:" << std::endl;
		std::cerr << error.what() << std::endl;
		std::cerr << "Interrupting code execution" << std::endl;
		std::abort();
	}
    
    double avg = lattice.get_average_number_of_adjacent_transpositions_per_bond();
    
    std::cout << "Number of bonds: " << lattice.get_nbonds() << std::endl;
    std::cout << "average number of adja transpos per bond: " << avg << std::endl;
	
	std::cout << "leaving main" << std::endl;
	
	return 0;
}
