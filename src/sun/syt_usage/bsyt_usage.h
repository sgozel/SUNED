// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_BSYT_USAGE_H
#define SUN_BSYT_USAGE_H

#include <vector>
#include <utility>

#include "../syt/bsyt.h"
#include "../irrep/irrep.h"
#include "../../common/datatypes.h"

namespace sun {

std::vector<tbSYT> get_SYT(const Irrep & alpha);

std::vector<int> get_column(const tbSYT syt, const Irrep& alpha);

std::pair<int, int> get_column_k_k_plus_one(const tbSYT syt, const unsigned int k);

int get_axial_distance(const tbSYT y, const std::vector<int> & cy, const int i, const int j);

void fullsimplify_development(std::vector<tbSYT> & ydev,
                              std::vector<double> & coeffdev,
                              const UINT64 sorted_up_to = 0);

void develop_consecutive_number(const Irrep & alpha,
								const std::vector<tbSYT> & ydev_in,
                                const std::vector<double> & coeffdev_in,
                                std::vector<tbSYT> & ydev_out,
                                std::vector<double> & coeffdev_out,
                                const int & k);

void develop_consecutive_number_inplace(const Irrep & alpha,
										std::vector<tbSYT> & ydev,
										std::vector<double> & coeffdev,
										const int & k);

} // namespace sun

#endif

