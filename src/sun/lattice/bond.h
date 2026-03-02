// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_BOND_H
#define SUN_BOND_H

#include <vector>
#include <string>
#include <cstdint>

#include "../permutation/permutation.h"

namespace sun {

struct Bond {
	Bond() {
		couplingValue = 0.0;
		couplingName = std::string("");
		
		id_ = 0;
		cacheit_ = false;
		parent_ = -1;
		child_ = -1;
	}
	
    double couplingValue;
    std::string couplingName;
    Permutation p;
    std::vector<AdjacentTransposition> ops;
    
    uint16_t id_;
    bool cacheit_;
    int parent_;
    int child_;
};

}

#endif
