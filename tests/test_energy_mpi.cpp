// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <string>
#include <mpi.h>

#include "nlohmann/json.hpp"

#include "testenergyutils.h"
#include "version.h"
#include "sun/heisenberg_mpi/hb_fund_matrix_engine_mpi.h"


TEST(HBFundMatrixEngineMPI, EnergySU3Chain)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	
	for (const auto& testsample : samples) {
		nlohmann::json inputParam = {
			{"N", testsample.N}, 
			{"Ns", testsample.Ns}, 
			{"alpha", testsample.alpha.get_vector()}, 
			{"latticefile", testsample.lattice}, 
			{"J", 1.0}, 
			{"max_iter", 100}, 
			{"tol_ritz", 1.0e-14}, 
			{"tol_residual", 1.0e-14}, 
			{"logging", false}, 
			{"checkpointing", false}, 
			{"dump_matrices", false}
		};
		
		sun::HBFundMatrixEngineMPI engine(inputParam);
		engine.initEngine();
		engine.build_matrix_lookups();
		double energy = engine.eig("multiply_mpi_matrix_v1");
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}


TEST(HBFundMatrixEngineMPI, NumaEnergySU3Chain)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	
	for (const auto& testsample : samples) {	
		nlohmann::json inputParam = {
			{"N", testsample.N}, 
			{"Ns", testsample.Ns}, 
			{"alpha", testsample.alpha.get_vector()}, 
			{"latticefile", testsample.lattice}, 
			{"J", 1.0}, 
			{"max_iter", 100}, 
			{"tol_ritz", 1.0e-14}, 
			{"tol_residual", 1.0e-14}, 
			{"logging", false}, 
			{"checkpointing", false}, 
			{"dump_matrices", false}
		};
		
		sun::HBFundMatrixEngineMPI engine(inputParam);
		engine.initEngine();
		engine.build_matrix_lookups();
		double energy = engine.eig("multiply_mpi_matrix_v1_numa");
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}


TEST(HBFundMatrixEngineMPI, NumaCheckpointing)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	EnergySample testsample = samples.back();

	nlohmann::json inputParam1 = {
			{"N", testsample.N}, 
			{"Ns", testsample.Ns}, 
			{"alpha", testsample.alpha.get_vector()}, 
			{"latticefile", testsample.lattice}, 
			{"J", 1.0}, 
			{"max_iter", 17}, 
			{"tol_ritz", 1.0e-14}, 
			{"tol_residual", 1.0e-14}, 
			{"logging", false}, 
			{"checkpointing", true}, 
			{"dump_matrices", false}
	};
	
	sun::HBFundMatrixEngineMPI engine1(inputParam1);
	engine1.initEngine();
	engine1.build_matrix_lookups();
	double energy1 = engine1.eig("multiply_mpi_matrix_v1_numa");
	
	nlohmann::json inputParam2 = {
			{"N", testsample.N}, 
			{"Ns", testsample.Ns}, 
			{"alpha", testsample.alpha.get_vector()}, 
			{"latticefile", testsample.lattice}, 
			{"J", 1.0}, 
			{"max_iter", 100}, 
			{"tol_ritz", 1.0e-14}, 
			{"tol_residual", 1.0e-14}, 
			{"logging", false}, 
			{"checkpointing", true}, 
			{"dump_matrices", false}
	};
	
	sun::HBFundMatrixEngineMPI engine2(inputParam2);
	engine2.initEngine();
	engine2.build_matrix_lookups();
	double energy2 = engine2.eig("multiply_mpi_matrix_v1_numa");
	
	EXPECT_NEAR(energy1, testsample.energy, 1e-1);
	EXPECT_NEAR(energy2, testsample.energy, 1e-12);
}
