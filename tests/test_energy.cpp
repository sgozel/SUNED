// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <string>

#include "nlohmann/json.hpp"

#include "testenergyutils.h"
#include "version.h"
#include "sun/heisenberg/hb_fund_mf_engine.h"
#include "sun/heisenberg/hb_fund_matrix_engine.h"


TEST(HBFundMatrixEngine, EnergySU3Chain)
{
	PRINT_SUNED_VERSION
	
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
		double energy = engine.eig("multiply_v1_openmp");
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}


TEST(HBFundMatrixFreeEngine, EnergySU3Chain)
{
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
		double energy = engine.eig("multiply_v1_openmp");
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}
