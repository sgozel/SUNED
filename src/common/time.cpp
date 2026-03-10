// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "time.h"

#include <iostream>
#include <iomanip>

#ifdef SG_USE_MPI
#include <mpi.h>
#endif

void time(
	const std::chrono::time_point<std::chrono::high_resolution_clock>& tstart, 
	const std::chrono::time_point<std::chrono::high_resolution_clock>& tend,
	const std::string& def)

{
	std::chrono::duration<double, std::milli> dt = tend - tstart;
    double elapsed = dt.count();
    std::ios_base::fmtflags coutflags(std::cout.flags());
#ifdef SG_USE_MPI
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	if (mpi_rank == 0) {
#endif
    std::cout << std::fixed;
    std::cout << std::setprecision(3);
    std::cout << def << " time = " << std::setw(9) << std::right 
			  << elapsed/1000 << " s" << std::endl;
    std::cout.flags(coutflags);
#ifdef SG_USE_MPI
	}
#endif
}


void time(
	const std::chrono::time_point<std::chrono::high_resolution_clock>& tstart, 
	const std::string& def)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> tend = std::chrono::high_resolution_clock::now();
	time(tstart, tend, def);
}
