// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef LANCZOSMPI_H
#define LANCZOSMPI_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include <mpi.h>

#include "../tmatrix/tmatrix.h"
#include "../lanczos/lanczosparams.h"

namespace lanczosmpi {


std::vector<double> residual(Tmatrix & tmat, const unsigned int k);


std::vector<double> ritz_value_stabilization(Tmatrix & tmat, const unsigned int k);


bool convergence(Tmatrix & tmat, const LanczosParams & lp);


void verify_convergence(Tmatrix & tmat, const unsigned int cpt, const LanczosParams & lp, const bool isConverged);


// y <--- y + a*x
template<class coeff_t>
void axpy(std::vector<coeff_t>& y, coeff_t a, const std::vector<coeff_t>& x);


// Perform one Lanczos step with 3 Lanczos vectors
template <class coeff_t, class type_mult>
inline void lanczos_step(std::vector<coeff_t> &, 
						 std::vector<coeff_t> &, 
						 std::vector<coeff_t> &,
                         double &, 
                         double &, 
                         type_mult);


// Lanczos diagonalization with 3 Lanczos vectors -- no output eigenvector
template <class coeff_t, class type_mult, class type_conv>
Tmatrix lanczos(const type_mult & multiply, 
				const type_conv & converge, 
				const UINT64 dimension, 
				const LanczosParams & lp);


// Lanczos diagonalization with 3 Lanczos vectors with output of eigenvector
template <class coeff_t, class type_mult, class type_conv>
Tmatrix lanczos_eigvec(const type_mult & multiply, 
					   const type_conv & converge, 
					   const UINT64 dimension, 
					   std::vector<coeff_t>& eigvec, 
					   const LanczosParams & lp);

#include "../lanczos/lanczos_utils.hpp"

#include "lanczos_checkpointing_mpi.hpp"

#include "lanczos_three_vectors_mpi.hpp"

}

#endif
