// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef MPI_UTILS_HPP
#define MPI_UTILS_HPP

#include <complex>
#include <mpi.h>

#include "datatypes.h"

template<class T> MPI_Datatype mpi_type();

template<> inline MPI_Datatype mpi_type<float>()  { return MPI_FLOAT; }
template<> inline MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }
template<> inline MPI_Datatype mpi_type<std::complex<float>>()  { return MPI_COMPLEX; }
template<> inline MPI_Datatype mpi_type<std::complex<double>>() { return MPI_DOUBLE_COMPLEX; }
template<> inline MPI_Datatype mpi_type<int16_t>()   { return MPI_INT16_T; }
template<> inline MPI_Datatype mpi_type<uint16_t>()   { return MPI_UINT16_T; }
template<> inline MPI_Datatype mpi_type<int32_t>()   { return MPI_INT32_T; }
template<> inline MPI_Datatype mpi_type<uint32_t>()   { return MPI_UINT32_T; }
template<> inline MPI_Datatype mpi_type<int64_t>()   { return MPI_INT64_T; }
template<> inline MPI_Datatype mpi_type<uint64_t>()   { return MPI_UINT64_T; }

#endif
