// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef LANCZOSPARAMS_H
#define LANCZOSPARAMS_H

#include <iostream>
#include <iomanip>
#include <string>

#ifdef SG_USE_MPI
#include <mpi.h>
#endif

#include "../nlohmann/json.hpp"

#ifdef SG_USE_MPI
namespace lanczosmpi {
#else
namespace lanczos {
#endif

struct LanczosParams {
	
	LanczosParams()
	{
#ifdef SG_USE_MPI
		int mpi_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
		
		// Lanczos init vector
		seed = 42;
#ifdef SG_USE_MPI
		seed += mpi_rank;
#endif
		
		// Lanczos iterations
		tol_residual = 1.0e-12;
		tol_ritz = 1.0e-12;
		min_iter = 1;
		max_iter = 1000;
		k = 1;
		conv_check_freq = 1;
		
		// Eigenvectors
		dump_eigvec = true;
		eigvec_folder = std::string("");
		
		// Logging
		logging = true;
		logging_full = false;
		logging_folder = std::string("");
		log_freq = 1;
		
		// Checkpointing
		checkpointing = false;
		checkpoint_folder = std::string("");
#ifdef SG_USE_MPI
		checkpoint_file = std::string("checkpoint_rank") + std::to_string(mpi_rank) + ".bin";
#else
		checkpoint_file = std::string("checkpoint.bin");
#endif
		checkpoint_frequency = 5;
		keep_checkpoint = false;
	}
	

	LanczosParams(nlohmann::json const& inputParam)
	{
		
#ifdef SG_USE_MPI
		int mpi_rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
#endif
		
		// Lanczos init vector
		seed = inputParam.value("seed", 42);
#ifdef SG_USE_MPI
		seed += mpi_rank;
#endif
		
		// Lanczos iterations
		tol_residual = inputParam.value("tol_residual", 1.0e-12);
		tol_ritz = inputParam.value("tol_ritz", 1.0e-12);
		min_iter = inputParam.value("min_iter", 1);
		max_iter = inputParam.value("max_iter", 1000);
		k = inputParam.value("k", 1);
		conv_check_freq = inputParam.value("convergence_check_frequency", 1);
		
		// Eigenvector
		dump_eigvec = inputParam.value("dump_eigvec", true);
		eigvec_folder = inputParam.value("eigvec_folder_path", std::string(""));
		
		// Logging
		logging = inputParam.value("logging", true);
		logging_full = inputParam.value("logging_full", false);
		logging_folder = inputParam.value("logging_folder_path", "");
		if ((!logging_folder.empty()) && (logging_folder.back()!='/')) {
			logging_folder += std::string("/");
		}
		log_freq = inputParam.value("logging_frequency", 1);
		
		// Checkpointing
		checkpointing = inputParam.value("checkpointing", false);
		checkpoint_folder = inputParam.value("checkpoint_folder_path", "");
		if ((!checkpoint_folder.empty()) && (checkpoint_folder.back()!='/')) {
			checkpoint_folder += std::string("/");
		}
		checkpoint_file = inputParam.value("checkpoint_file", "checkpoint.bin");
#ifdef SG_USE_MPI
		auto dot_pos = checkpoint_file.rfind('.');
		if (dot_pos != std::string::npos) {
			checkpoint_file = checkpoint_file.substr(0, dot_pos) 
							+ "_rank" + std::to_string(mpi_rank) 
							+ checkpoint_file.substr(dot_pos);
		} else {
			checkpoint_file = checkpoint_file + "_rank" + std::to_string(mpi_rank);
		}
#endif
		checkpoint_file = checkpoint_folder + checkpoint_file;
		checkpoint_frequency = inputParam.value("checkpoint_frequency", 5);
		keep_checkpoint = inputParam.value("keep_checkpoint", false);
	}
	
	void print() {
		std::ios_base::fmtflags coutflags(std::cout.flags());
		std::cout << "------------------------------------------" << std::endl;
		std::cout << "Lanczos parameters:" << std::endl;
		std::cout << "------------------------------------------" << std::endl;
		std::cout << "seed ------------------- : " << seed << std::endl;
		//std::cout << "min_iter --------------- : " << min_iter << std::endl;
		std::cout << "max_iter --------------- : " << max_iter << std::endl;
		std::cout << "tol_residual ----------- : " << std::scientific << std::setprecision(3) << tol_residual << std::endl;
		std::cout << "tol_ritz --------------- : " << std::scientific << std::setprecision(3) << tol_ritz << std::endl;
		//std::cout << "k ---------------------- : " << k << std::endl;
		std::cout << "conv_check_freq -------- : " << conv_check_freq << std::endl;
#ifdef SG_LANCZOS_EIGVEC
		std::cout << "dump_eigvec ------------ : " << dump_eigvec << std::endl;
		if (dump_eigvec == true) {
			std::cout << "eigvec_folder ---------- : " << eigvec_folder << std::endl;
		}
#endif
		std::cout << "logging ---------------- : " << logging << std::endl;
		if (logging == true) {
			std::cout << "logging_folder --------- : " << logging_folder << std::endl;
			std::cout << "log_freq --------------- : " << log_freq << std::endl;
		}
		std::cout << "checkpointing ---------- : " << checkpointing << std::endl;
		if (checkpointing == true) {
			std::cout << "checkpoint_folder ------ : " << checkpoint_folder << std::endl;
			std::cout << "checkpoint_file -------- : " << checkpoint_file << std::endl;
			std::cout << "checkpoint_frequency --- : " << checkpoint_frequency << std::endl;
			std::cout << "keep_checkpoint -------- : " << keep_checkpoint << std::endl;
		}
		std::cout.flags(coutflags);
	}
	
	// Initialization of Lanczos vector
	unsigned int seed;
	
	// Convergence/Iteration parameters
    unsigned int min_iter;
    unsigned int max_iter;
    double tol_residual;
    double tol_ritz;
    unsigned int k;
    unsigned int conv_check_freq;
    
    // eigenvector
    bool dump_eigvec;
    std::string eigvec_folder;
    
    // Logging parameters
    bool logging;
	bool logging_full;
    unsigned int log_freq;
    std::string logging_folder;
    
    // Checkpointing parameters
    bool checkpointing;
    std::string checkpoint_folder;
    std::string checkpoint_file;
    unsigned int checkpoint_frequency;
    bool keep_checkpoint;
};

} // namespace lanczos(mpi)

#endif
