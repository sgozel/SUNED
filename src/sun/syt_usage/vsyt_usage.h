// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_VSYT_USAGE_H
#define SUN_VSYT_USAGE_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <limits>
#include <cmath>
#include <stdexcept>

#include "../syt/vsyt.h"
#include "../irrep/irrep.h"
#include "../utils/young_factor.h"
#include "../utils/utils.h"
#include "../../common/datatypes.h"


namespace sun {

template<class type_t>
std::vector<vSYT<type_t>> get_SYT(const Irrep & alpha, const UINT64 from = 0, UINT64 dim = std::numeric_limits<UINT64>::max());


template<class type_t>
std::vector<int> get_column(const vSYT<type_t>& syt, const Irrep& alpha);


template<class type_t>
std::pair<int, int> get_column_k_k_plus_one(const vSYT<type_t> syt, const unsigned int k);


template<class type_t>
int get_axial_distance(const vSYT<type_t> & y, const std::vector<int> & cy, const int i, const int j);


template<class type_t>
void fullsimplify_development(std::vector<vSYT<type_t>> & ydev,
                              std::vector<double> & coeffdev,
                              const UINT64 sorted_up_to = 0);


template<class type_t>
void develop_consecutive_number(const Irrep & alpha,
								const std::vector<vSYT<type_t>> & ydev_in,
                                const std::vector<double> & coeffdev_in,
                                std::vector<vSYT<type_t>> & ydev_out,
                                std::vector<double> & coeffdev_out,
                                const int & k);


template<class type_t>
void develop_consecutive_number_inplace(const Irrep & alpha,
										std::vector<vSYT<type_t>> & ydev,
										std::vector<double> & coeffdev,
										const int & k);

#include "vsyt_usage.hpp"

} // namespace sun

#endif
