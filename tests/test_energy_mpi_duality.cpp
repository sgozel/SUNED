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


TEST(SYSTEMDUALITY, EnergySU5Chain)
{
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
		PRINT_SUNED_VERSION
	}
	
	const std::string mvm_method("multiply_mpi_matrix_v1");
	
	std::string filename("TEST_DATA_ENERGY_SU5.json");
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
		engine.init();
		engine.build_matrix_lookups();
		double energy = engine.eigenvalue(mvm_method);
		
		EXPECT_NEAR(energy, testsample.energy, 1e-12);
	}
}
