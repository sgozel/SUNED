// Copyright 2026 Samuel GOZEL, GNU GPLv3
#ifndef MPI_COMM_HPP
#define MPI_COMM_HPP

#include <iostream>
#include <vector>
#include <mpi.h>

#include "mpi_utils.hpp"



template<class coeff_t>
/*void alltoallv_chunked(const std::vector<coeff_t>& send,
					   const std::vector<int64_t>& sendcounts,
					   const std::vector<int64_t>& sdispls,
                       std::vector<coeff_t>& recv,
                       const std::vector<int64_t>& recvcounts,
                       const std::vector<int64_t>& rdispls,
                       MPI_Comm comm = MPI_COMM_WORLD)*/
void alltoallv_chunked(const coeff_t* send,
                       const std::vector<int64_t>& sendcounts,
                       const std::vector<int64_t>& sdispls,
                       coeff_t* recv,
                       const std::vector<int64_t>& recvcounts,
                       const std::vector<int64_t>& rdispls,
                       MPI_Comm comm = MPI_COMM_WORLD)
{	
    int world_size;
    MPI_Comm_size(comm, &world_size);

	/*
	// NOT POSSIBLE ANYMORE, BECAUSE recv IS A POINTER
    // Verify recv buffer is large enough
    int64_t total_recvcounts = 0;
    for (int r = 0; r < world_size; ++r) {
        total_recvcounts += recvcounts[r];
    }
    if (static_cast<int64_t>(recv.size()) < total_recvcounts) {
        std::cerr << "alltoallv_chunked: recv buffer too small (size="
                  << recv.size() << ", required=" << total_recvcounts << ")"
                  << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }*/

    const int64_t chunk_size = std::numeric_limits<int>::max();
    
    int64_t total_send = sdispls[world_size-1] + sendcounts[world_size-1];
	int64_t total_recv = rdispls[world_size-1] + recvcounts[world_size-1];
    
    int one_chunk = ((total_send <= chunk_size) && (total_recv <= chunk_size)) ? 1 : 0;
	int all_one_chunk;
	MPI_Allreduce(&one_chunk, &all_one_chunk, 1, MPI_INT, MPI_MIN, comm);
    
    if (all_one_chunk) {
		// single chunk - no reallocation needed - only cast counts/displs to int
		std::vector<int> send_counts(world_size);
		std::vector<int> recv_counts(world_size);
        std::vector<int> send_displs(world_size);
        std::vector<int> recv_displs(world_size);
        for (int r = 0; r < world_size; ++r) {
            send_counts[r] = static_cast<int>(sendcounts[r]);
            recv_counts[r] = static_cast<int>(recvcounts[r]);
            send_displs[r] = static_cast<int>(sdispls[r]);
            recv_displs[r] = static_cast<int>(rdispls[r]);
        }
		
        MPI_Alltoallv(
            send, //send.data(),
            send_counts.data(),
            send_displs.data(),
            mpi_type<coeff_t>(),
            recv, //recv.data(),
            recv_counts.data(),
            recv_displs.data(),
            mpi_type<coeff_t>(),
            comm
        );
        
	} else {
		// Several chunks are needed
		
		std::cerr << "PROBLEM: alltoallv : several chunks needed. Need a review" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
		
		int64_t max_send = *std::max_element(sendcounts.begin(), sendcounts.end());
		int64_t max_recv = *std::max_element(recvcounts.begin(), recvcounts.end());
		int64_t max_count = std::max(max_send, max_recv);
		
		int64_t n_chunks = (max_count + chunk_size - 1) / chunk_size;
		MPI_Allreduce(MPI_IN_PLACE, &n_chunks, 1, MPI_INT64_T, MPI_MAX, comm);
		
		std::vector<int> send_counts(world_size, 0);
		std::vector<int> recv_counts(world_size, 0);
		std::vector<int> send_displs(world_size, 0);
		std::vector<int> recv_displs(world_size, 0);
		
		// Pack send staging buffer
		std::vector<coeff_t> send_buf(std::numeric_limits<int>::max());
		std::vector<coeff_t> recv_buf(std::numeric_limits<int>::max());

		for (int64_t chunk_id = 0; chunk_id < n_chunks; ++chunk_id)
		{
			const int64_t offset = chunk_id * chunk_size;
			
			int64_t total_send = 0;
			int64_t total_recv = 0;

			for (int r = 0; r < world_size; ++r)
			{
				const int64_t rem_send = std::max(int64_t(0), sendcounts[r] - offset);
				const int64_t rem_recv = std::max(int64_t(0), recvcounts[r] - offset);
				send_counts[r] = static_cast<int>(std::min(rem_send, chunk_size));
				recv_counts[r] = static_cast<int>(std::min(rem_recv, chunk_size));
				send_displs[r] = static_cast<int>(total_send);
				recv_displs[r] = static_cast<int>(total_recv);
				total_send += send_counts[r];
				total_recv += recv_counts[r];
			}
			
			for (int r = 0; r < world_size; ++r) {
				/*
				// vector version
				std::copy(send.begin() + sdispls[r] + offset,
						  send.begin() + sdispls[r] + offset + send_counts[r],
						  send_buf.begin() + send_displs[r]);*/
				// pointer version
				std::copy(send + sdispls[r] + offset,
						  send + sdispls[r] + offset + send_counts[r],
						  send_buf.begin() + send_displs[r]);
			}
			
			MPI_Alltoallv(
				send_buf.data(),
				send_counts.data(),
				send_displs.data(),
				mpi_type<coeff_t>(),
				recv_buf.data(),
				recv_counts.data(),
				recv_displs.data(),
				mpi_type<coeff_t>(),
				comm
			);

			// Unpack recv staging buffer
			for (int r = 0; r < world_size; ++r) {
				/*
				// vector version
				std::copy(recv_buf.begin() + recv_displs[r],
						  recv_buf.begin() + recv_displs[r] + recv_counts[r],
						  recv.begin() + rdispls[r] + offset);*/
				// Pointer version
				std::copy(recv_buf.begin() + recv_displs[r],
						  recv_buf.begin() + recv_displs[r] + recv_counts[r],
						  recv + rdispls[r] + offset);
			}
		}
	}
}


template<class coeff_t>
void alltoallv(/*const std::vector<coeff_t>& send,
               const std::vector<int64_t>& sendcounts,
               const std::vector<int64_t>& sdispls,
               std::vector<coeff_t>& recv,
               const std::vector<int64_t>& recvcounts,
               const std::vector<int64_t>& rdispls,
               MPI_Comm comm = MPI_COMM_WORLD)*/
			   const coeff_t* send,
               const std::vector<int64_t>& sendcounts,
               const std::vector<int64_t>& sdispls,
               coeff_t* recv,
               const std::vector<int64_t>& recvcounts,
               const std::vector<int64_t>& rdispls,
               MPI_Comm comm = MPI_COMM_WORLD)
{
#if MPI_VERSION >= 4
    static_assert(sizeof(MPI_Count) == sizeof(int64_t), "MPI_Count size mismatch");

    MPI_Alltoallv_c(
        send,
        sendcounts.data(),
        sdispls.data(),
        mpi_type<coeff_t>(),
        recv,
        recvcounts.data(),
        rdispls.data(),
        mpi_type<coeff_t>(),
        comm
    );
#else
    alltoallv_chunked(
		send,
		sendcounts,
		sdispls,
		recv,
		recvcounts,
		rdispls,
		comm);
#endif
}


// thin wrapper - takes vectors, delegates to pointer version
template<class coeff_t>
void alltoallv(const std::vector<coeff_t>& send,
               const std::vector<int64_t>& sendcounts,
               const std::vector<int64_t>& sdispls,
               std::vector<coeff_t>& recv,
               const std::vector<int64_t>& recvcounts,
               const std::vector<int64_t>& rdispls,
               MPI_Comm comm = MPI_COMM_WORLD)
{
    alltoallv(
		send.data(),
		sendcounts,
		sdispls,
		recv.data(),
		recvcounts,
		rdispls,
		comm
	);
}

#endif
