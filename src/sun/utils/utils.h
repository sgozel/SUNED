// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_UTILS_H
#define SUN_UTILS_H

#include <vector>

#include "../../common/datatypes.h"
#include "../irrep/irrep.h"

namespace sun {

void reduce_by_divide(std::vector<unsigned int>& numvec, std::vector<unsigned int>& denomvec);

UINT64 multiplicity(const Irrep & alpha);

template<class SYT>
void print_syt(const SYT& syt, const Irrep & alpha, const bool & zeroBased=true);

#include "utils.hpp"

} // namespace sun

#endif
