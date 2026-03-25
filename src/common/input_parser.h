// Copyright 2026 Samuel GOZEL, GNU GPLv3
#ifndef SUNED_INPUT_PARSER_H
#define SUNED_INPUT_PARSER_H

#include <iostream>
#include <fstream>
#include <string>

#ifdef SG_USE_MPI
#include <mpi.h>
#endif

#include "nlohmann/json.hpp"
using json = nlohmann::json;


json parseInputArguments(int argc, char* argv[])
{
#ifdef SG_USE_MPI
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
	
	std::ifstream inputFileStream;
	json inputParam;
	
	try {
		if (argc == 1) {
#ifdef SG_USE_MPI
			if (mpi_rank == 0) {
				std::cerr << "Missing input parameter file / parameters. Aborting." << std::endl;
			}
			MPI_Abort(MPI_COMM_WORLD, 1);
#else
			std::cerr << "Missing input parameter file / parameters. Aborting." << std::endl;
			std::abort();
#endif
		} else {
			int st = 1;
			std::string firstArg = argv[st];
			if (firstArg.substr(0, 2) != "--") {
				// Attemp to read from a .json file				
				inputFileStream.open(firstArg);
				if (inputFileStream.is_open()) {
					try {
						inputFileStream >> inputParam;
						inputFileStream.close();
						st += 1;
					} catch (...) {
#ifdef SG_USE_MPI
						if (mpi_rank == 0) {
							std::cerr << "Failed to parse JSON file: " << firstArg << ". Aborting." << std::endl;
						}
						MPI_Abort(MPI_COMM_WORLD, 1);
#else
						std::cerr << "Failed to parse JSON file: " << firstArg << ". Aborting." << std::endl;
						std::abort();
#endif
					}
				} else {
					// File not found - typo in file name / file path
#ifdef SG_USE_MPI
					if (mpi_rank == 0) {
						std::cerr << "Failed to open JSON file: " << firstArg << ". Aborting." << std::endl;
					}
					MPI_Abort(MPI_COMM_WORLD, 1);
#else
					std::cerr << "Failed to open JSON file: " << firstArg << ". Aborting." << std::endl;
					std::abort();
#endif
				}
			}
			
			// Catch a dangling flag with no value
			if ((argc-st) % 2 == 1) {
#ifdef SG_USE_MPI
				if (mpi_rank == 0) {
					std::cerr << "Wrong parity of input arguments. Aborting." << std::endl;
				}
				MPI_Abort(MPI_COMM_WORLD, 1);
#else
				std::cerr << "Wrong parity of input arguments. Aborting." << std::endl;
				std::abort();
#endif
			}
			
			// Parse optional --key value pairs
			for ( ; st < argc - 1; st += 2) {
				std::string flag = argv[st];
				if (flag.substr(0, 2) == "--") {
					
					std::string key = flag.substr(2);
					std::string value = argv[st + 1];
					
					if (value.front() == '[') {
						inputParam[key] = json::parse(value);
					} else {
						try {
							size_t pos;
							int ivalue = std::stoi(value, &pos);
							if (pos == value.size()) {
								inputParam[key] = ivalue;  // fully consumed --> integer
							} else {
								inputParam[key] = std::stod(value); // partially consumed --> double
							}
						} catch (...) {
							if (value == "true") {
								inputParam[key] = true;
							} else if (value == "false") {
								inputParam[key] = false;
							} else {
								// string
								inputParam[key] = value;
							}
						}
					}
				} else {
#ifdef SG_USE_MPI
					if (mpi_rank == 0) {
						std::cerr << "Unexpected argument: " << flag << ". Expected --key value format." << std::endl;
					}
					MPI_Abort(MPI_COMM_WORLD, 1);
#else
					std::cerr << "Unexpected argument: " << flag << ". Expected --key value format." << std::endl;
					std::abort();
#endif
				}
			}
		}
	}
	catch (...) {
#ifdef SG_USE_MPI
		if (mpi_rank == 0) {
			std::cerr << "Caught exception at top level in [main]." << std::endl;
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
#else
		std::cerr << "Caught exception at top level in [main]." << std::endl;
		std::abort();	
#endif
	}
	
	return inputParam;
}

#endif
