// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef LANCZOSMPI_H
#define LANCZOSMPI_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <omp.h>
#include <mpi.h>

#include "../common/numa.h"
#include "../tmatrix/tmatrix.h"
#include "../lanczos/lanczosparams.h"

namespace lanczosmpi {

template<class coeff_t>
void lanczos_init_vector(sg_vec<coeff_t>& v, const UINT64 dimension, const unsigned int seed);


template<class coeff_t>
void numa_lanczos_init_vector(sg_vec<coeff_t>& v, const UINT64 dimension, const unsigned int seed);


std::vector<double> residual(Tmatrix & tmat, unsigned int k);


std::vector<double> ritz_value_stabilization(Tmatrix & tmat, unsigned int k);


bool convergence(Tmatrix & tmat, const LanczosParams & lp);


void verify_convergence(Tmatrix & tmat, const unsigned int cpt, const LanczosParams & lp, const bool isConverged);


template<typename E, typename coeff_t, class Alloc>
void dump_eigpair(const E energy,
				  const std::vector<coeff_t, Alloc>& eigvec,
				  const unsigned int index, 
				  const LanczosParams& lp);


template<typename E, typename coeff_t, class Alloc>
bool load_eigpair(std::pair<E, std::vector<coeff_t, Alloc>>& eigpair, 
				  const unsigned int index, 
				  const LanczosParams& lp);


// y <--- y + a*x
template<class coeff_t>
void axpy(sg_vec<coeff_t>& y, coeff_t a, const sg_vec<coeff_t>& x);


// NUMA-Aware y <--- y + a*x
template<class coeff_t>
void numa_axpy(sg_vec<coeff_t>& y, coeff_t a, const sg_vec<coeff_t>& x);


// Perform one Lanczos step with 3 Lanczos vectors
template <class coeff_t, class type_mult>
inline void lanczos_step(sg_vec<coeff_t> &, 
						 sg_vec<coeff_t> &, 
						 sg_vec<coeff_t> &,
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
					   sg_vec<coeff_t>& eigvec, 
					   const LanczosParams & lp);

#include "../lanczos/lanczos_utils.hpp"

#include "../lanczos/lanczos_checkpointing.hpp"

#include "lanczos_three_vectors_mpi.hpp"

}

#endif
