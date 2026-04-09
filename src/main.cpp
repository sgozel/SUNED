// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>
#include <utility>

#include "nlohmann/json.hpp"
using json = nlohmann::json;

#include "version.h"
#include "common/input_parser.h"
#include "common/numa.h"
#include "sun/irrep/irrep.h"
#include "sun/heisenberg/hb_fund_matrix_engine.h"


int main(int argc, char* argv[])
{
	PRINT_SUNED_VERSION
	
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    // READ INPUT PARAMETERS
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    
    json inputParam = parseInputArguments(argc, argv);
	
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	// RUN SIMULATION
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

	sun::HBFundMatrixEngine engine(inputParam);
    engine.init();
    engine.build_matrix_lookups();
    
    const std::string mvm_method("multiply_v1_openmp");

    #ifdef SG_EIGVALS
    // Compute eigenvalue
    double energy = engine.eigenvalue(mvm_method);
    std::cout << "In main: from eigenvalue(), energy = " << energy << std::endl;
    #endif

    #ifdef SG_EIGVECS
    // Compute eigenpair
    std::pair<double, sg_vec<double>> eigpair = engine.eigenpair(mvm_method);
    std::cout << "In main: from eigenpair(), energy = " << eigpair.first << std::endl;
    engine.check_eigvec(eigpair, mvm_method);
	#endif
	//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	
	std::cout << "leaving main" << std::endl;
	
	return 0;
}
