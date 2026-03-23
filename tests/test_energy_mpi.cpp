// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <mpi.h>

#include "nlohmann/json.hpp"

#include "testenergyutils.h"
#include "testutils.h"
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
	for (size_t i=0; i < 1; ++i) {
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


TEST(HBFundMatrixEngineMPI, EigenpairSU3Chain)
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
	for (size_t i=0; i < 1; ++i) {
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


TEST(HBFundMatrixEngineMPI, MatrixDumpingLoading)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	const std::string mvm_method("multiply_mpi_matrix_v1");
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	EnergySample testsample = samples[5];

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
			{"dump_matrices", true},
			{"matrix_dump_folder_path", "."}
	};
	
	// engine which builds the matrix lookups
	sun::HBFundMatrixEngineMPI engine1(inputParam);
	engine1.init();
	engine1.build_matrix_lookups();
	
	// engine which reads the matrix lookups and computes energy
	sun::HBFundMatrixEngineMPI engine2(inputParam);
	engine2.init();
	engine2.build_matrix_lookups();
	double energy = engine2.eigenvalue(mvm_method);
	
	EXPECT_NEAR(energy, testsample.energy, 1e-12);
}


TEST(HBFundMatrixEngineMPI, Checkpointing)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	const std::string mvm_method("multiply_mpi_matrix_v1");
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	EnergySample testsample = samples[5];

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


TEST(HBFundMatrixEngineMPI, DumpLoadEigvec)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	const std::string mvm_method("multiply_mpi_matrix_v1");
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	EnergySample testsample = samples[5];

	nlohmann::json inputParam = {
			{"N", testsample.N}, 
			{"Ns", testsample.Ns}, 
			{"alpha", testsample.alpha.get_vector()}, 
			{"latticefile", testsample.lattice}, 
			{"J", 1.0}, 
			{"logging", false}, 
			{"checkpointing", false}, 
			{"dump_matrices", false},
			{"dump_eigvec", true}
	};
	
	// Compute eigenpair and dump it to disk
	sun::HBFundMatrixEngineMPI engine1(inputParam);
	engine1.init();
	engine1.build_matrix_lookups();
	std::pair<double, sg_vec<double>> eigpair_computed = engine1.eigenpair(mvm_method);
	
	// Load eigenpair from disk
	sun::HBFundMatrixEngineMPI engine2(inputParam);
	std::pair<double, sg_vec<double>> eigpair_loaded = engine2.load_eigpair(0);
	
	EXPECT_EQ(eigpair_computed.first, eigpair_loaded.first);
	EXPECT_EQ(eigpair_computed.second.size(), eigpair_loaded.second.size());
	expect_equal_vec(eigpair_computed.second, eigpair_loaded.second);
}


TEST(HBFundMatrixEngineMPI, Correlations)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	const std::string mvm_method("multiply_mpi_matrix_v1");
	
	std::string filename("TEST_DATA_ENERGY_SU3.json");
	std::vector<EnergySample> samples = read_energy_test_data(filename);
	EnergySample testsample = samples[5];

	nlohmann::json inputParam = {
			{"N", testsample.N}, 
			{"Ns", testsample.Ns}, 
			{"alpha", testsample.alpha.get_vector()}, 
			{"latticefile", testsample.lattice}, 
			{"J", 1.0}, 
			{"tol_ritz", 1.0e-14}, 
			{"tol_residual", 1.0e-14}, 
			{"logging", false}, 
			{"checkpointing", false}, 
			{"dump_matrices", false},
			{"dump_eigvec", false}
	};
	
	// Compute eigenpair
	sun::HBFundMatrixEngineMPI engine1(inputParam);
	engine1.init();
	engine1.build_matrix_lookups();
	std::pair<double, sg_vec<double>> eigpair = engine1.eigenpair(mvm_method);
	
	// Check eigenvalue
	EXPECT_NEAR(eigpair.first, testsample.energy, 1e-12);
	
	// Compute correlations
	std::vector<double> C = engine1.correlations(eigpair.second);
	
	EXPECT_EQ(C[0], 1.0); // quite a trivial test at this point ...
}
