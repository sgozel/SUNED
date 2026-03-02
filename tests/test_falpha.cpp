// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <gtest/gtest.h>
#include <fstream>
#include <vector>
#include <stdexcept>

#include "nlohmann/json.hpp"

#include "version.h"
#include "common/datatypes.h"
#include "sun/irrep/irrep.h"
#include "sun/utils/utils.h"


struct MultiplicitySample {
    sun::Irrep alpha;
    UINT64 falpha;
};

std::vector<MultiplicitySample> read_multiplicity_test_data()
{
	std::ifstream file("TEST_DATA_MULTIPLICITY.json");
    if (!file.is_open()) {
        throw std::runtime_error("Error: could not open file.");
    }
    
    nlohmann::json data;
    file >> data;
    
    std::vector<MultiplicitySample> falpha_samples;
    for (const auto& item : data) {
        MultiplicitySample e;
        std::vector<unsigned int> temp_alpha = item["alpha"].get<std::vector<unsigned int>>();
        e.alpha = sun::Irrep(temp_alpha);
        e.falpha = item["falpha"].get<UINT64>();
        falpha_samples.push_back(e);
    }
    return falpha_samples;
}

TEST(SUNUtils, Multiplicity)
{
	PRINT_SUNED_VERSION
	
	std::vector<MultiplicitySample> samples = read_multiplicity_test_data();
	
	for (const auto& testsample : samples)
	{	
		UINT64 falpha = sun::multiplicity(testsample.alpha);
		EXPECT_EQ(falpha, testsample.falpha);
	}
}
