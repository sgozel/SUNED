// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
#include "sun/lattice/lattice.h"


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
	
	if (!inputParam.contains("latticefile")) {
		throw std::runtime_error("Missing latticefile in input .json file.");
	}
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
