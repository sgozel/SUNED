// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <string>
#include <utility>

#include "nlohmann/json.hpp"

#include "testenergyutils.h"
#include "version.h"
#include "common/numa.h"
#include "sun/heisenberg/hb_fund_mf_engine.h"
#include "sun/heisenberg/hb_fund_matrix_engine.h"


TEST(HBFundMatrixEngine, EnergySU3Chain)
{
	PRINT_SUNED_VERSION
	
	const std::string mvm_method("multiply_v1_openmp");
	
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
		
		sun::HBFundMatrixEngine engine(inputParam);
		engine.initEngine();
		engine.build_matrix_lookups();
		double energy = engine.eigenvalue(mvm_method);
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}


TEST(HBFundMatrixEngine, NumaEnergySU3Chain)
{
	PRINT_SUNED_VERSION
	
	const std::string mvm_method("multiply_v1_openmp_numa");
	
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
		
		sun::HBFundMatrixEngine engine(inputParam);
		engine.initEngine();
		engine.build_matrix_lookups();
		double energy = engine.eigenvalue(mvm_method);
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}



TEST(HBFundMatrixEngine, NumaEigenpairSU3Chain)
{
	PRINT_SUNED_VERSION
	
	const std::string mvm_method("multiply_v1_openmp_numa");
	
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
		
		sun::HBFundMatrixEngine engine(inputParam);
		engine.initEngine();
		engine.build_matrix_lookups();
		std::pair<double, sg_vec<double>> eigpair = engine.eigenpair(mvm_method);
		
		double energy = eigpair.first;
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
		
		double distance = engine.check_eigvec(eigpair, mvm_method);
		EXPECT_LT(distance, 1.0e-6);
	}
}



TEST(HBFundMatrixEngine, NumaCheckpointing)
{
	PRINT_SUNED_VERSION
	
	const std::string mvm_method("multiply_v1_openmp_numa");
	
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
	
	sun::HBFundMatrixEngine engine1(inputParam1);
	engine1.initEngine();
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
	
	sun::HBFundMatrixEngine engine2(inputParam2);
	engine2.initEngine();
	engine2.build_matrix_lookups();
	double energy2 = engine2.eigenvalue(mvm_method);
	
	EXPECT_NEAR(energy1, testsample.energy, 1e-1);
	EXPECT_NEAR(energy2, testsample.energy, 1e-12);
}


TEST(HBFundMatrixFreeEngine, EnergySU3Chain)
{
	PRINT_SUNED_VERSION
	
	const std::string mvm_method("multiply_v1_openmp");
	
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
			{"checkpointing", false}
		};
		
		sun::HBFundMatrixFreeEngine engine(inputParam);
		engine.initEngine();
		double energy = engine.eigenvalue(mvm_method);
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}
