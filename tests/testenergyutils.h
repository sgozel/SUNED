// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <iostream>
#include <vector>
#include <string>
#ifdef SG_USE_MPI
#include <mpi.h>
#else
#include <stdexcept>
#endif

#include "nlohmann/json.hpp"

#include "sun/irrep/irrep.h"


struct EnergySample {
    unsigned int N;
    unsigned int Ns;
    sun::Irrep alpha;
    std::string lattice;
    std::string BC;
    double energy;
};

std::vector<EnergySample> read_energy_test_data(const std::string& filename)
{
	std::ifstream file(filename);
    if (!file.is_open()) {
#ifdef SG_USE_MPI        
        std::cerr << "Error: could not open test data file " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
#else
		throw std::runtime_error("Error: could not open test data file " + filename);
#endif
    }
    
    nlohmann::json data;
    file >> data;
    
    std::vector<EnergySample> energy_samples;
    for (const auto& item : data) {
        EnergySample e;
        e.N = item["N"].get<unsigned int>();
        e.Ns = item["Ns"].get<unsigned int>();
        std::vector<unsigned int> temp_alpha = item["alpha"].get<std::vector<unsigned int>>();
        e.alpha = sun::Irrep(temp_alpha);
        e.lattice = item["lattice"].get<std::string>();
        e.BC = item["BC"].get<std::string>();
        e.lattice = std::string("HB_") + e.lattice + "_" + std::to_string(e.Ns) + '_' + e.BC + ".lattice";
        e.energy = item["energy"].get<double>();
        energy_samples.push_back(e);
    }
    return energy_samples;
}

#endif

