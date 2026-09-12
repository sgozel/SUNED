// Copyright 2026 Samuel GOZEL, GNU GPLv3
#ifndef MPI_COMM_HPP
#define MPI_COMM_HPP

#include <iostream>
#include <vector>
#include <mpi.h>

#include "mpi_utils.hpp"

namespace mpi {

namespace detail {


// MPI 4 Isend/Irecv communication using large counts MPI_Isend_c / MPI_Irecv_c
template<class coeff_t>
void isend_irecv_large(const coeff_t* send,
                       const std::vector<int64_t>& sendcounts,
                       const std::vector<int64_t>& sdispls,
                       coeff_t* recv,
                       const std::vector<int64_t>& recvcounts,
                       const std::vector<int64_t>& rdispls,
                       MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);

    std::vector<MPI_Request> reqs;
    reqs.reserve(2 * world_size);  // upper bound

    // Post all receives
    for (int rank = 0; rank < world_size; ++rank) {
        if (recvcounts[rank] > 0) {
            MPI_Request req;
            MPI_Irecv_c(recv + rdispls[rank],  // receive buffer
                        static_cast<MPI_Count>(recvcounts[rank]), // count
                        mpi_type<coeff_t>(),   // datatype
                        rank,                  // source process
						0,                     // tag
						comm, 
						&req);
            reqs.push_back(req);
        }
    }

    // Post all sends
    for (int rank = 0; rank < world_size; ++rank) {
        if (sendcounts[rank] > 0) {
            MPI_Request req;
            MPI_Isend_c(send + sdispls[rank],   // send buffer
                        static_cast<MPI_Count>(sendcounts[rank]), // count
                        mpi_type<coeff_t>(),    // datatype
                        rank,  					// destination process
						0, 						// tag
						comm, 
						&req);
            reqs.push_back(req);
        }
    }

    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
}


// MPI 3 implementation using chunked Isend/Irecv communication
template<class coeff_t>
void isend_irecv_chunked(const coeff_t* send,
                         const std::vector<int64_t>& sendcounts,
                         const std::vector<int64_t>& sdispls,
                         coeff_t* recv,
                         const std::vector<int64_t>& recvcounts,
                         const std::vector<int64_t>& rdispls,
                         MPI_Comm comm)
{
    int world_size;
    MPI_Comm_size(comm, &world_size);

    // Max elements per chunk — stay well below INT_MAX
    static constexpr int64_t chunk_size = 1 << 30;  // ~1 billion elements ~ 8.6GB for doubles

    std::vector<MPI_Request> reqs;
    reqs.reserve(2 * world_size);  // will grow if chunking kicks in

    // Post all receives
    for (int rank = 0; rank < world_size; ++rank) {
        if (recvcounts[rank] > 0) {
            int64_t remaining = recvcounts[rank];
            int64_t offset    = rdispls[rank];
            int     tag       = 0;
            while (remaining > 0) {
                const int count = static_cast<int>(std::min(remaining, chunk_size));
                MPI_Request req;
                MPI_Irecv(recv + offset,
                          count,
                          mpi_type<coeff_t>(),
                          rank, 
						  tag, 
						  comm, 
						  &req);
                reqs.push_back(req);
                offset    += count;
                remaining -= count;
                tag += 1;
            }
        }
    }

    // Post all sends
    for (int rank = 0; rank < world_size; ++rank) {
        if (sendcounts[rank] > 0) {
            int64_t remaining = sendcounts[rank];
            int64_t offset    = sdispls[rank];
            int     tag       = 0;
            while (remaining > 0) {
                const int count = static_cast<int>(std::min(remaining, chunk_size));
                MPI_Request req;
                MPI_Isend(send + offset,
                          count,
                          mpi_type<coeff_t>(),
                          rank, 
						  tag, 
						  comm, 
						  &req);
                reqs.push_back(req);
                offset    += count;
                remaining -= count;
                tag += 1;
            }
        }
    }

    MPI_Waitall(static_cast<int>(reqs.size()), reqs.data(), MPI_STATUSES_IGNORE);
}


template<class coeff_t>
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

    const int64_t chunk_size = 1 << 30; // std::numeric_limits<int>::max();
    
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
            send,
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
		// Several chunks are needed - rely on isend_irecv_chunked
		isend_irecv_chunked(send, sendcounts, sdispls,
                           	recv, recvcounts, rdispls, comm);
	}
}


} // namespace detail


// public API - takes send and receive buffers as pointers
template<class coeff_t>
void alltoallv(const coeff_t* send,
               const std::vector<int64_t>& sendcounts,
               const std::vector<int64_t>& sdispls,
               coeff_t* recv,
               const std::vector<int64_t>& recvcounts,
               const std::vector<int64_t>& rdispls,
               MPI_Comm comm = MPI_COMM_WORLD)
{
#if MPI_VERSION >= 4
    static_assert(sizeof(MPI_Count) == sizeof(int64_t), "MPI_Count size mismatch");
    static_assert(sizeof(MPI_Aint)  == sizeof(int64_t), "MPI_Aint size mismatch");

    int world_size;
    MPI_Comm_size(comm, &world_size);

    std::vector<MPI_Count> sc(world_size);
    std::vector<MPI_Count> rc(world_size);
    std::vector<MPI_Aint>  sd(world_size);
    std::vector<MPI_Aint>  rd(world_size);
    for (int r = 0; r < world_size; ++r) {
        sc[r] = static_cast<MPI_Count>(sendcounts[r]);
        rc[r] = static_cast<MPI_Count>(recvcounts[r]);
        sd[r] = static_cast<MPI_Aint>(sdispls[r]);
        rd[r] = static_cast<MPI_Aint>(rdispls[r]);
    }

    MPI_Alltoallv_c(
        send,
        sc.data(),
        sd.data(),
        mpi_type<coeff_t>(),
        recv,
        rc.data(),
        rd.data(),
        mpi_type<coeff_t>(),
        comm
    );
#else
    detail::alltoallv_chunked(
		send,
		sendcounts,
		sdispls,
		recv,
		recvcounts,
		rdispls,
		comm);
#endif
}


// public API, thin wrapper - takes send and receive buffers as vectors, delegates to pointer version
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


// public API - equivalent to alltoallv, implemented with Isend/Irecv
template<class coeff_t>
void isend_irecv(const coeff_t* send,
                 const std::vector<int64_t>& sendcounts,
                 const std::vector<int64_t>& sdispls,
                 coeff_t* recv,
                 const std::vector<int64_t>& recvcounts,
                 const std::vector<int64_t>& rdispls,
                 MPI_Comm comm = MPI_COMM_WORLD)
{
#if MPI_VERSION >= 4
    static_assert(sizeof(MPI_Count) == sizeof(int64_t), "MPI_Count size mismatch");
    detail::isend_irecv_large(send, sendcounts, sdispls,
                              recv, recvcounts, rdispls, comm);
#else
    detail::isend_irecv_chunked(send, sendcounts, sdispls,
                                recv, recvcounts, rdispls, comm);
#endif
}

} // namespace mpi

#endif
