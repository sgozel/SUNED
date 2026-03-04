// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef LANCZOS_H
#define LANCZOS_H

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

#include "../tmatrix/tmatrix.h"
#include "lanczosparams.h"

namespace lanczos {

template<class coeff_t>
void lanczos_init_vector(std::vector<coeff_t>& v, const UINT64 dimension, const unsigned int seed);


std::vector<double> residual(Tmatrix & tmat, const unsigned int k);


std::vector<double> ritz_value_stabilization(Tmatrix & tmat, const unsigned int k);


bool convergence(Tmatrix & tmat, const LanczosParams & lp);


void verify_convergence(Tmatrix & tmat, const unsigned int cpt, const LanczosParams & lp, const bool isConverged);


template<typename coeff_t>
void dump_eigvec(const std::vector<coeff_t>& eigvec, const unsigned int index, const LanczosParams & lp);


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


// Perform one Lanczos step with 2 Lanczos vectors
template <class coeff_t, class type_mult>
inline void lanczos_step_two_vectors(std::vector<coeff_t> &, 
						 std::vector<coeff_t> &, 
                         double &, 
                         double &, 
                         type_mult);


// Lanczos diagonalization with 2 Lanczos vectors -- no output eigenvector
template <class coeff_t, class type_mult, class type_conv>
Tmatrix lanczos_two_vectors(const type_mult & multiply, 
							const type_conv & converge, 
							const UINT64 dimension, 
							const LanczosParams & lp);


// Lanczos diagonalization with 2 Lanczos vectors with output of eigenvector
template <class coeff_t, class type_mult, class type_conv>
Tmatrix lanczos_two_vectors_eigvec(const type_mult & multiply, 
								   const type_conv & converge, 
								   const UINT64 dimension, 
								   std::vector<coeff_t> & eigvec, 
								   const LanczosParams & lp);


#include "lanczos_utils.hpp"

#include "lanczos_checkpointing.hpp"

#include "lanczos_two_vectors.hpp"

#include "lanczos_three_vectors.hpp"

} // namespace lanczos

#endif
