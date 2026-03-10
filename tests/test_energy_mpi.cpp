// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <mpi.h>

#include "nlohmann/json.hpp"

#include "testenergyutils.h"
#include "version.h"
#include "common/numa.h"
#include "sun/heisenberg_mpi/hb_fund_matrix_engine_mpi.h"


TEST(HBFundMatrixEngineMPI, EnergySU3Chain)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	const std::string mvm_method("multiply_mpi_matrix_v1");
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	
	//for (const auto& testsample : samples) {
	for (size_t i=0; i < 4; ++i) {
		const auto & testsample = samples[i];
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
		engine.init();
		engine.build_matrix_lookups();
		double energy = engine.eigenvalue(mvm_method);
		
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
	
	const std::string mvm_method("multiply_mpi_matrix_v1_numa");
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	
	//for (const auto& testsample : samples) {
	for (size_t i=0; i < 4; ++i) {
		const auto & testsample = samples[i];
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
		engine.init();
		engine.build_matrix_lookups();
		double energy = engine.eigenvalue(mvm_method);
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}



TEST(HBFundMatrixEngineMPI, NumaEigenpairSU3Chain)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	const std::string mvm_method("multiply_mpi_matrix_v1_numa");
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	
	//for (const auto& testsample : samples) {
	for (size_t i=0; i < 4; ++i) {
		const auto & testsample = samples[i];
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
		engine.init();
		engine.build_matrix_lookups();
		std::pair<double, sg_vec<double>> eigpair = engine.eigenpair(mvm_method);
		
		double energy = eigpair.first;
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
		
		double distance = engine.check_eigvec(eigpair, mvm_method);
		EXPECT_LT(distance, 1.0e-6);
	}
}


TEST(HBFundMatrixEngineMPI, NumaCheckpointing)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	const std::string mvm_method("multiply_mpi_matrix_v1_numa");
	
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
	engine1.init();
	engine1.build_matrix_lookups();
	double energy1 = engine1.eigenvalue(mvm_method);
	
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
	engine2.init();
	engine2.build_matrix_lookups();
	double energy2 = engine2.eigenvalue(mvm_method);
	
	EXPECT_NEAR(energy1, testsample.energy, 1e-1);
	EXPECT_NEAR(energy2, testsample.energy, 1e-12);
}
