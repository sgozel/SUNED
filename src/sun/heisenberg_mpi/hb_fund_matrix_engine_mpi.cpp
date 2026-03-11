// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "hb_fund_matrix_engine_mpi.h"

#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <chrono>
#include <omp.h>
#include <mpi.h>

#include "../utils/utils.h"
#include "../../common/time.h"
#include "../../common/mpi_utils.hpp"
#include "../../common/mpi_comm.hpp"

#ifdef SG_USE_BASIC_SYT
#include "../syt_usage/vsyt_usage.h"
#else
#include "../syt_usage/bsyt_usage.h"
#endif


namespace sun {

HBFundMatrixEngineMPI::HBFundMatrixEngineMPI(nlohmann::json const& inputParam)
: HBFundEngineMPI(inputParam)
{	
	dump_matrices_ = inputParam.value("dump_matrices", false);
	
	if (dump_matrices_ == true) {
		if (!inputParam.contains("matrix_dump_folder_path")) {
			std::cerr << "Missing matrix_dump_folder_path in input .json file for HBFundEngineMPI." << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		matrix_dump_path_ = inputParam["matrix_dump_folder_path"];
		if (matrix_dump_path_.back()!='/') {
			matrix_dump_path_ += std::string("/");
		}
	}
	
	P_.resize(alpha_.n()-1);
	mpi_offdiag_nodes_.resize(alpha_.n()-1);
	mpi_offdiag_nodes_acc_.resize(alpha_.n()-1);
	for (unsigned int k=0; k<alpha_.n()-1; ++k) {
		mpi_offdiag_nodes_[k].resize(mpi_world_size_);
		std::fill(mpi_offdiag_nodes_[k].begin(), mpi_offdiag_nodes_[k].end(), 0);
		mpi_offdiag_nodes_acc_[k].resize(mpi_world_size_);
		std::fill(mpi_offdiag_nodes_acc_[k].begin(), mpi_offdiag_nodes_acc_[k].end(), 0);
	}
	mpi_nb_offdiag_.resize(alpha_.n()-1);
	std::fill(mpi_nb_offdiag_.begin(), mpi_nb_offdiag_.end(), 0);
	
	mpi_local_index_base_.resize(alpha_.n()-1);
	mpi_local_index_friend_.resize(alpha_.n()-1);
	
	Y_bounds_.resize(mpi_world_size_);
	
	if (mpi_rank_==0) {
		std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << "N = " << N_ << std::endl;
		std::cout << "Ns = " << Ns_ << std::endl;
		std::cout << "Target irrep: " << std::endl;
		alpha_.print();
		std::cout << std::endl;
		lanczosparams_.print();
		lattice_.print_sites();
		lattice_.print_bonds();
	}
}


void HBFundMatrixEngineMPI::init()
{
	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	
	dimension_ = multiplicity(alpha_);
	
	if (mpi_rank_==0) {
		std::cout << "dimension = " << dimension_ << std::endl;
	}
	
	mpi_get_local_dimension();
	print_mpi_details();
	
	// Each MPI process only generates its own SYTs
	
	const UINT64 from = mpi_rank_ * mpi_bare_dimension_;
	
	#ifdef SG_USE_BASIC_SYT
	std::cerr << "Currently, BASIC_SYT is unsupported on MPI application." << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
	Y_ = get_SYT<SYTel>(alpha_, from, mpi_dimension_);
	#else
	Y_ = get_SYT(alpha_, from, mpi_dimension_);
	#endif
	
	// Each process sends its Y_[0] to all processes
	const SYT_value_t first_val = Y_[0].value();
	std::vector<SYT_value_t> all_first(mpi_world_size_);
	MPI_Allgather(
		&first_val,              // send buffer
		1,                       // send count
		mpi_type<SYT_value_t>(), // send type
		all_first.data(),        // recv buffer
		1,                       // recv count per process
		mpi_type<SYT_value_t>(), // recv type
		MPI_COMM_WORLD
	);

	// Each process sends its Y_[Y_.size()-1] to all processes
	const SYT_value_t last_val  = Y_[Y_.size()-1].value();
	std::vector<SYT_value_t> all_last(mpi_world_size_);
	MPI_Allgather(
		&last_val,               // send buffer
		1,                       // send count
		mpi_type<SYT_value_t>(), // send type
		all_last.data(),         // recv buffer
		1,                       // recv count per process
		mpi_type<SYT_value_t>(), // recv type
		MPI_COMM_WORLD
	);
	
	for (int rank = 0; rank < mpi_world_size_; ++rank) {
		Y_bounds_[rank].first = SYT(all_first[rank]);
		Y_bounds_[rank].second = SYT(all_last[rank]);
	}
	
	work_.resize(mpi_dimension_);
#ifdef SG_USE_NUMA
	#pragma omp parallel for schedule(static)
	for (UINT64 i = 0; i < mpi_dimension_; ++i) {
		work_[i] = 0.0;
	}
#endif
	
	if (mpi_rank_==0) {
		double factor = 1e6;
		std::string units("MB");
		if (8*Y_.size()>=1e9) {
			factor *= 1000;
			units = std::string("GB");
		}
		
		#ifdef SG_USE_BASIC_SYT
		double memAllSYTs = alpha_.n() * sizeof(SYTel) * static_cast<double>(dimension_)/factor;
		double memLocalSYTs = alpha_.n() * sizeof(SYTel) * static_cast<double>(mpi_dimension_)/factor;
		#else
		double memAllSYTs = sizeof(SYT_value_t) * static_cast<double>(dimension_)/factor;
		double memLocalSYTs = sizeof(SYT_value_t) * static_cast<double>(mpi_dimension_)/factor;
		#endif
		double memMatrixLookups = sizeof(typePk) * (alpha_.n()-1) * static_cast<double>(mpi_dimension_)/factor;
		double memLanczos = 3 * sizeof(double) * static_cast<double>(mpi_dimension_)/factor;
		double memWorkArray = sizeof(double) * static_cast<double>(mpi_dimension_)/factor;
		double memTotal = memMatrixLookups + memLanczos + memWorkArray;
		
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << ":::::::: TOTAL MEMORY REQUIREMENTS ::::::::" << std::endl;
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << "Total SYTs Memory: " << memAllSYTs << units << std::endl;
		std::cout << "Local SYTs Memory: " << memLocalSYTs << units << std::endl;
		std::cout << "Matrix lookups:    " << memMatrixLookups << units << std::endl;
		std::cout << "3 Lanczos vectors: " << memLanczos << units << std::endl;
		std::cout << "work_ array:       " << memWorkArray << units << std::endl;
		std::cout << "--------------" << std::endl;
		std::cout << "Total: " << memTotal << units << std::endl;
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
	}
	
	time(t0, std::string("init"));
}


void HBFundMatrixEngineMPI::build_matrix_lookups()
{	
	std::chrono::time_point<std::chrono::high_resolution_clock> t_start = std::chrono::high_resolution_clock::now();
	
	// one more security measure
	static_assert(std::is_trivially_copyable<SYT>::value, 
		"SYT must be trivially copyable for MPI byte transfer");
	
	std::vector<UINT64> offdiag_indices(mpi_dimension_);
	
	for (unsigned int k=0; k<alpha_.n()-1; ++k) {
		
		std::chrono::time_point<std::chrono::high_resolution_clock> tk_start = std::chrono::high_resolution_clock::now();
		std::string transpo_string = std::string("(") + std::to_string(k) + ", " + std::to_string(k+1) + ")";
		
		std::fill(offdiag_indices.begin(), offdiag_indices.end(), 0);
		
		P_[k].resize(mpi_dimension_);
		
		// Attempt to load P_[k] from file if it exists
		if ((dump_matrices_==true) && (load_matrix(k))) {
			if (mpi_rank_ == 0) {
				time(tk_start, std::string("Load from file ")+transpo_string);
			}
            continue;
        }
		
		// For each remote rank, a list of local_i indices and yfriend SYTs
		std::vector<std::vector<UINT64>> pending_i(mpi_world_size_);
		std::vector<std::vector<SYT_value_t>> pending_syt(mpi_world_size_);
		
		#pragma omp parallel for schedule(guided)
		for (UINT64 i = 0; i < mpi_dimension_; ++i)
		{	
			const int rowk = Y_[i].get(k);
			const int rowkk = Y_[i].get(k+1);
			
			if (rowk == rowkk) {
				P_[k][i] = 1;
			} else {
				const std::pair<int, int> cy = get_column_k_k_plus_one(Y_[i], k);
				if (cy.first == cy.second) {
					P_[k][i] = -1;
				} else {		
					SYT yfriend = Y_[i];
					yfriend.exchange(k, k+1);
					
					// Extract the rank of the friend SYT
					int rankfriend;
					for (int rank = 0; rank < mpi_world_size_; ++rank) {
						if ((yfriend>=Y_bounds_[rank].first) && (yfriend<=Y_bounds_[rank].second)) {
							rankfriend = rank;
							break;
						}
					}
					
					if (rankfriend == mpi_rank_) {
						// The friend SYT belongs to this rank - we can search within the local collection of SYTs
						
						auto it = std::lower_bound(Y_.begin(), Y_.end(), yfriend);
						UINT64 index = it - Y_.begin(); // local index
						index += mpi_rank_ * mpi_bare_dimension_; // global index
						
						offdiag_indices[i] = index + 1; // add 1 to differentiate from value 0 used for diagonal elements
					
						// search to which process <index> belongs to
						const unsigned int rankfriend_v2 = mpi_rank_from_index(index);
						
						if (rankfriend_v2 != rankfriend) {
							std::cerr << "PROBLEM: missmatch between two different methods of extracting the rank of the friend SYT." << std::endl;
							MPI_Abort(MPI_COMM_WORLD, 1);
						}
						
						#pragma omp atomic
						mpi_offdiag_nodes_[k][rankfriend] += 1;

						#pragma omp atomic
						mpi_nb_offdiag_[k] += 1;
						
						// axial distance from k to k+1 in SYT Y_[i]
						// count +1 for each step made downwards or to the left
						// count -1 for each step made upwards or to the right
						const typePk ax = cy.first - rowk - cy.second + rowkk;
						
						P_[k][i] = -ax;
						
					} else {
						// register a request for a index extraction on a friend rank
						#pragma omp critical
						{
							pending_i[rankfriend].emplace_back(i);
							pending_syt[rankfriend].emplace_back(yfriend.value());
						}
					}
				}
			}
		} // end for i (local Hilbert space)
		
		// MPI EXCHANGE for remote friends
		
		// How many requests this rank is sending to each other rank
		std::vector<int64_t> send_counts(mpi_world_size_);
		for (int r = 0; r < mpi_world_size_; ++r) {
			send_counts[r] = pending_i[r].size();
		}

		// How many requests this rank will receive from each other rank
		std::vector<int64_t> recv_counts(mpi_world_size_);
		MPI_Alltoall(
			send_counts.data(),
			1,
			mpi_type<int64_t>(),
			recv_counts.data(),
			1,
			mpi_type<int64_t>(),
			MPI_COMM_WORLD
		);
		
		// Compute send displacements
		std::vector<int64_t> send_displs(mpi_world_size_, 0);
		for (int r = 1; r < mpi_world_size_; ++r) {
			send_displs[r] = send_displs[r-1] + send_counts[r-1];
		}
		const int64_t total_send = send_displs[mpi_world_size_-1] + send_counts[mpi_world_size_-1];
		
		// Compute recv displacements
		std::vector<int64_t> recv_displs(mpi_world_size_, 0);
		for (int r = 1; r < mpi_world_size_; ++r) {
			recv_displs[r] = recv_displs[r-1] + recv_counts[r-1];
		}
		const int64_t total_recv = recv_displs[mpi_world_size_-1] + recv_counts[mpi_world_size_-1];
		
		// Pack data to be sent into a contiguous buffers
		std::vector<SYT_value_t> send_buf_syt(total_send);
		for (int r = 0; r < mpi_world_size_; ++r) {
			std::copy(pending_syt[r].begin(),
					  pending_syt[r].end(),
					  send_buf_syt.begin() + send_displs[r]);
		}
		
		// free memory of pending_syt
		{ std::vector<std::vector<SYT_value_t>>().swap(pending_syt); }
		
		// receive buffer
		std::vector<SYT_value_t> recv_buf_syt(total_recv);
		
		// Exchange incoming friend SYT
		alltoallv(
			send_buf_syt,
			send_counts,
			send_displs,
			recv_buf_syt,
			recv_counts,
			recv_displs,
			MPI_COMM_WORLD
		);

		// search received SYTs on this local collection
		std::vector<UINT64> reply_buf(total_recv);

		#pragma omp parallel for schedule(guided)
		for (int64_t j = 0; j < total_recv; ++j) {
			const SYT yfriend(recv_buf_syt[j]);
			auto it = std::lower_bound(Y_.begin(), Y_.end(), yfriend);
			const UINT64 local_index = it - Y_.begin();
			const UINT64 global_index = mpi_rank_ * mpi_bare_dimension_ + local_index;
			reply_buf[j] = global_index;
		}
		
		// Send back the results of the requests to the original processes
		std::vector<UINT64> result_buf(total_send); // answers coming back to this process
		
		alltoallv(
			reply_buf,
			recv_counts,
			recv_displs,
			result_buf,
			send_counts,
			send_displs,
			MPI_COMM_WORLD
		);
		
		// FILL IN PASS
		
		for (int r = 0; r < mpi_world_size_; ++r)
		{
			#pragma omp parallel for schedule(static)
			for (int64_t j = 0; j < send_counts[r]; ++j)
			{
				const UINT64 local_i      = pending_i[r][j];
				const UINT64 global_index = result_buf[send_displs[r] + j];

				offdiag_indices[local_i] = global_index + 1; // +1 convention, same as local case

				const SYT& yi = Y_[local_i];
				const std::pair<int, int> cy = get_column_k_k_plus_one(yi, k);
				const int rowk  = yi.get(k);
				const int rowkk = yi.get(k+1);
				const typePk ax = cy.first - rowk - cy.second + rowkk;
				P_[k][local_i] = -ax;

				const unsigned int rankfriend = mpi_rank_from_index(global_index);
				#pragma omp atomic
				mpi_offdiag_nodes_[k][rankfriend] += 1;
				#pragma omp atomic
				mpi_nb_offdiag_[k] += 1;
			}
		}
		
		// mpi_nb_offdiag_[k] = number of off-diagonal element for transposition (k, k+1) found by this process
		// 
		// mpi_offdiag_nodes_[k][rank] = for transposition (k, k+1), number of images belonging to <rank>
		// 
		// offdiag_indices[i] = 0 for diagonal elements; index+1 of image SYT (k, k+1) Y_[i] (index is the index, from 0 to dimension_)
		
		for (int rank=0; rank<mpi_world_size_-1; ++rank) {
			mpi_offdiag_nodes_acc_[k][rank+1] = mpi_offdiag_nodes_acc_[k][rank] + mpi_offdiag_nodes_[k][rank];
		}
		
		UINT64 total_offdiags_count = mpi_offdiag_nodes_acc_[k][mpi_world_size_-1] + mpi_offdiag_nodes_[k][mpi_world_size_-1];
		
		if (total_offdiags_count != mpi_nb_offdiag_[k]) {
			std::cerr << "PROBLEM A" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		// array which will keep track of the number of friend SYTs belonging to each rank
		std::vector<UINT64> friend_nodes_count(mpi_world_size_, 0);
		
		// array which will keep track of the local index from this process for all off-diagonal elements
		mpi_local_index_base_[k].resize(total_offdiags_count);
		
		// array which will keep track of the local index of the image SYT in the friend process, for all off-diagonal elements
		std::vector<UINT64> mpi_local_index_friend_temp(total_offdiags_count, 0);
		
		UINT64 cpt = 0;
		
		for (UINT64 i=0; i<mpi_dimension_; ++i) {
			if (offdiag_indices[i]>0) {
				// extract index of image SYT
				const UINT64 index = offdiag_indices[i] - 1; // global index
				
				// compute the MPI rank which owns this SYT
				const unsigned int rankfriend = mpi_rank_from_index(index);
				
				// local index in this process of the original (base) SYT
				mpi_local_index_base_[k][mpi_offdiag_nodes_acc_[k][rankfriend] + friend_nodes_count[rankfriend]] = i;
				
				// local index of image (friend) SYT in friend process
				mpi_local_index_friend_temp[mpi_offdiag_nodes_acc_[k][rankfriend] + friend_nodes_count[rankfriend]] = mpi_local_index_from_global_index(index);
				
				friend_nodes_count[rankfriend] += 1;
				cpt += 1;
			}
		}
		
		UINT64 temp_test = 0;
		for (int rank=0; rank<mpi_world_size_; ++rank) {
			temp_test += friend_nodes_count[rank];
		}
		
		if (cpt!=total_offdiags_count) {
			std::cerr << "PROBLEM B" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		if (temp_test!=total_offdiags_count) {
			std::cerr << "PROBLEM C" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		std::vector<UINT64> friend_nodes_count_received(mpi_world_size_, 0);
		
		MPI_Alltoall(
			friend_nodes_count.data(),
			1,                                  // sendcount: 1 element per destination
			MPI_UINT64_T,
			friend_nodes_count_received.data(),
			1,                                  // recvcount: 1 element per source
			MPI_UINT64_T,
			MPI_COMM_WORLD
		);
		
		// check that friend_nodes_count[r] == friend_nodes_count_received[r]
		// meaning that each rank sends and receives as many elements from each other rank
		for (int r = 0; r < mpi_world_size_; ++r) {
			if (friend_nodes_count[r] != friend_nodes_count_received[r]) {
				std::cerr << "PROBLEM D: friend_nodes_count mismatch between rank " << mpi_rank_ 
						  << " and rank " << r
						  << " friend_nodes_count[" << r << "] = " << friend_nodes_count[r]
						  << "; friend_nodes_count_received[" << r << "] = " << friend_nodes_count_received[r]
						  << std::endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		
		mpi_local_index_friend_[k].resize(total_offdiags_count);
		
		// Because of the structure of the interaction, recvcounts will be the
		// same as sendcounts, and rdispls will be the same as sdispls
		
		const std::vector<int64_t>& sendrecvcounts = mpi_offdiag_nodes_[k];
		const std::vector<int64_t>& srdispls = mpi_offdiag_nodes_acc_[k];
		
		alltoallv(
			mpi_local_index_friend_temp,
			sendrecvcounts,
			srdispls,
			mpi_local_index_friend_[k], // receive buffer
			sendrecvcounts,
			srdispls,
			MPI_COMM_WORLD
		);
		
		time(tk_start, std::string("Building ")+transpo_string);
		
		if (dump_matrices_ == true) {
			dump_matrix(k);
		}
	}
	
	time(t_start, "build_matrix_lookups");
	
	free_basis();
}


void HBFundMatrixEngineMPI::free_basis()
{
	std::vector<SYT>().swap(Y_);
}


void HBFundMatrixEngineMPI::dump_matrix(const unsigned int k) const
{
    const std::string filename = matrix_dump_path_ 
							   + std::string("Pk_mpi_ws") + std::to_string(mpi_world_size_) 
							   + "_rank" + std::to_string(mpi_rank_) 
							   + "_k" + std::to_string(k) + ".bin";
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Cannot open file : " << filename << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Write mpi_dimension_
    out.write(reinterpret_cast<const char*>(&mpi_dimension_), sizeof(UINT64));

    // Write P_[k]
    out.write(reinterpret_cast<const char*>(P_[k].data()), mpi_dimension_ * sizeof(typePk));

    // Write mpi_nb_offdiag_[k]
    out.write(reinterpret_cast<const char*>(&mpi_nb_offdiag_[k]), sizeof(UINT64));

    // Write mpi_offdiag_nodes_[k]
    out.write(reinterpret_cast<const char*>(mpi_offdiag_nodes_[k].data()), mpi_world_size_ * sizeof(INT64));

    // Write mpi_offdiag_nodes_acc_[k]
    out.write(reinterpret_cast<const char*>(mpi_offdiag_nodes_acc_[k].data()), mpi_world_size_ * sizeof(INT64));

    // Write mpi_local_index_base_[k]
    const UINT64 base_size = mpi_local_index_base_[k].size();
    out.write(reinterpret_cast<const char*>(&base_size), sizeof(UINT64));
    out.write(reinterpret_cast<const char*>(mpi_local_index_base_[k].data()), base_size * sizeof(UINT64));

    // Write mpi_local_index_friend_[k]
    const UINT64 friend_size = mpi_local_index_friend_[k].size();
    out.write(reinterpret_cast<const char*>(&friend_size), sizeof(UINT64));
    out.write(reinterpret_cast<const char*>(mpi_local_index_friend_[k].data()), friend_size * sizeof(UINT64));
}


bool HBFundMatrixEngineMPI::load_matrix(const unsigned int k)
{
    const std::string filename = matrix_dump_path_ 
							   + std::string("Pk_mpi_ws") + std::to_string(mpi_world_size_) 
							   + "_rank" + std::to_string(mpi_rank_) 
							   + "_k" + std::to_string(k) + ".bin";
	
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        return false;
    }

    // Read and verify mpi_dimension_
    UINT64 stored_dimension = 0;
    in.read(reinterpret_cast<char*>(&stored_dimension), sizeof(UINT64));
    if (!in || stored_dimension != mpi_dimension_) {
        std::cerr << "Warning: dimension mismatch or read error in file: " << filename
                  << " (stored=" << stored_dimension << ", expected=" << mpi_dimension_ << ")" << std::endl;
        return false;
    }

    // Read P_[k]
    P_[k].resize(mpi_dimension_);
    in.read(reinterpret_cast<char*>(P_[k].data()), mpi_dimension_ * sizeof(typePk));
    if (!in) {
        std::cerr << "Warning: failed to read P_[k] from file: " << filename << std::endl;
        return false;
    }

    // Read mpi_nb_offdiag_[k]
    in.read(reinterpret_cast<char*>(&mpi_nb_offdiag_[k]), sizeof(UINT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_nb_offdiag_[k] from file: " << filename << std::endl;
        return false;
    }

    // Read mpi_offdiag_nodes_[k]
    mpi_offdiag_nodes_[k].resize(mpi_world_size_);
    in.read(reinterpret_cast<char*>(mpi_offdiag_nodes_[k].data()), mpi_world_size_ * sizeof(INT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_offdiag_nodes_[k] from file: " << filename << std::endl;
        return false;
    }

    // Read mpi_offdiag_nodes_acc_[k]
    mpi_offdiag_nodes_acc_[k].resize(mpi_world_size_);
    in.read(reinterpret_cast<char*>(mpi_offdiag_nodes_acc_[k].data()), mpi_world_size_ * sizeof(INT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_offdiag_nodes_acc_[k] from file: " << filename << std::endl;
        return false;
    }

    // Read mpi_local_index_base_[k]
    UINT64 base_size = 0;
    in.read(reinterpret_cast<char*>(&base_size), sizeof(UINT64));
    if (base_size != mpi_nb_offdiag_[k]) {
        std::cerr << "Warning: base_size mismatch in file: " << filename << std::endl;
        return false;
    }
    mpi_local_index_base_[k].resize(base_size);
    in.read(reinterpret_cast<char*>(mpi_local_index_base_[k].data()), base_size * sizeof(UINT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_local_index_base_[k] from file: " << filename << std::endl;
        return false;
    }

    // Read mpi_local_index_friend_[k]
    UINT64 friend_size = 0;
    in.read(reinterpret_cast<char*>(&friend_size), sizeof(UINT64));
    if (friend_size != mpi_nb_offdiag_[k]) {
        std::cerr << "Warning: friend_size mismatch in file: " << filename << std::endl;
        return false;
    }
    mpi_local_index_friend_[k].resize(friend_size);
    in.read(reinterpret_cast<char*>(mpi_local_index_friend_[k].data()), friend_size * sizeof(UINT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_local_index_friend_[k] from file: " << filename << std::endl;
        return false;
    }

    return true;
}


void HBFundMatrixEngineMPI::multiply(const sg_vec<double> & w, sg_vec<double> & u, const double & a, const std::string & method) const
{
    if (method=="multiply_mpi_matrix_v1") {
        multiply_mpi_matrix_v1(w, u, a);
	} else if (method=="multiply_mpi_matrix_v1_numa") {
        multiply_mpi_matrix_v1_numa(w, u, a);
	} else {
        std::cerr << "Multiply method undefined" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}


template <class coeff_t>
void HBFundMatrixEngineMPI::multiply_mpi_matrix_v1(const sg_vec<coeff_t>& w, sg_vec<coeff_t>& u, const double a) const
{
	// u <--- H*w - a*u
	
	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	
	std::for_each(u.begin(), u.end(), [a](coeff_t& el) {el*=(-a);});

	std::vector<int64_t> sendrecvcounts(mpi_world_size_, 0);
	std::vector<int64_t> srdispls(mpi_world_size_, 0);

	// one allocation before loop on all bonds and all transpositions in 
	// each bond. Allocation size based on the max number of off-diagonal
	// elements accross all transpositions
	const UINT64 max_offdiag = *std::max_element(mpi_nb_offdiag_.begin(), mpi_nb_offdiag_.end());
	std::vector<coeff_t> send_coeffs(max_offdiag);
	std::vector<coeff_t> recv_coeffs(max_offdiag);

	unsigned int cpt_bond = 0;

	for (const auto& bond : lattice_.bonds)
	{
		//std::chrono::time_point<std::chrono::high_resolution_clock> tbond_0 = std::chrono::high_resolution_clock::now();
		
		std::copy(w.begin(), w.end(), work_.begin());
		
		for (unsigned int j=0; j<bond.ops.size(); ++j)
		{
			//std::chrono::time_point<std::chrono::high_resolution_clock> tbondj_0 = std::chrono::high_resolution_clock::now();
			
			const unsigned int k = bond.ops[j].getk();
			
			for (int r=0; r<mpi_world_size_; ++r) {
				sendrecvcounts[r] = static_cast<int64_t>(mpi_offdiag_nodes_[k][r]);
				srdispls[r] = static_cast<int64_t>(mpi_offdiag_nodes_acc_[k][r]);
			}
			
			// gather all coefficients from this process which will be 
			// sent to all processes
			#pragma omp parallel for schedule(guided)
			for(UINT64 i=0; i<mpi_nb_offdiag_[k]; ++i) {
				const double rho = 1.0/static_cast<double>(P_[k][mpi_local_index_base_[k][i]]);
				send_coeffs[i] = work_[mpi_local_index_base_[k][i]] * std::sqrt(1.0-rho*rho);
			}
			
			//======================================
			// MPI communication of coefficients
			//======================================
			
			alltoallv(
				send_coeffs,
				sendrecvcounts,
				srdispls,
				recv_coeffs,
				sendrecvcounts,
				srdispls,
				MPI_COMM_WORLD
			);
			
			//======================================
			// Perform update of <work> array
			//======================================
			
			// All diagonal terms
			#pragma omp parallel for schedule(guided)
			for (UINT64 i=0; i<mpi_dimension_; ++i) {
				work_[i] *= 1.0/static_cast<double>(P_[k][i]);
			}
			
			// All off-diagonal terms - no risk of data race at this point
			#pragma omp parallel for schedule(guided)
			for (UINT64 i=0; i<mpi_nb_offdiag_[k]; ++i) {
				work_[mpi_local_index_friend_[k][i]] += recv_coeffs[i];
			}
		} // for j (operations in a bond)
		
		//===========================
		// UPDATE OF LANCZOS VECTOR
		//===========================
		
		const double J = bond.couplingValue;
		for (UINT64 i=0; i<mpi_dimension_; ++i) {
			u[i] += J * work_[i];
		}
		cpt_bond += 1;
	}
	
	time(t0, "multiply");
}


template <class coeff_t>
void HBFundMatrixEngineMPI::multiply_mpi_matrix_v1_numa(const sg_vec<coeff_t>& w, sg_vec<coeff_t>& u, const double a) const
{
	// u <--- H*w - a*u
	
	std::chrono::time_point<std::chrono::high_resolution_clock> t0 = std::chrono::high_resolution_clock::now();
	
	#pragma omp parallel for schedule(static)
	for (UINT64 i = 0; i < mpi_dimension_; ++i) {
		u[i] *= -a;
	}

	std::vector<int64_t> sendrecvcounts(mpi_world_size_, 0);
	std::vector<int64_t> srdispls(mpi_world_size_, 0);

	// one allocation before loop on all bonds and all transpositions in 
	// each bond. Allocation size based on the max number of off-diagonal
	// elements accross all transpositions
	const UINT64 max_offdiag = *std::max_element(mpi_nb_offdiag_.begin(), mpi_nb_offdiag_.end());
	
	std::vector<coeff_t> send_coeffs(max_offdiag); // use sg_vec instead
	std::vector<coeff_t> recv_coeffs(max_offdiag); // use sg_vec instead

	for (const auto& bond : lattice_.bonds)
	{
		#pragma omp parallel for schedule(static)
		for (UINT64 i = 0; i < mpi_dimension_; ++i) {
			work_[i] = w[i];
		}
		
		for (unsigned int j=0; j<bond.ops.size(); ++j)
		{
			const unsigned int k = bond.ops[j].getk();
			
			for (int r=0; r<mpi_world_size_; ++r) {
				sendrecvcounts[r] = static_cast<int64_t>(mpi_offdiag_nodes_[k][r]);
				srdispls[r] = static_cast<int64_t>(mpi_offdiag_nodes_acc_[k][r]);
			}
			
			// gather all coefficients from this process which will be 
			// sent to all processes
			//#pragma omp parallel for schedule(guided)
			#pragma omp parallel for schedule(static)
			for(UINT64 i=0; i<mpi_nb_offdiag_[k]; ++i) {
				const double rho = 1.0/static_cast<double>(P_[k][mpi_local_index_base_[k][i]]);
				send_coeffs[i] = work_[mpi_local_index_base_[k][i]] * std::sqrt(1.0-rho*rho);
			}
			
			//======================================
			// MPI communication of coefficients
			//======================================
			
			alltoallv(
				send_coeffs,
				sendrecvcounts,
				srdispls,
				recv_coeffs,
				sendrecvcounts,
				srdispls,
				MPI_COMM_WORLD
			);
			
			//======================================
			// Perform update of <work> array
			//======================================
			
			// All diagonal terms
			//#pragma omp parallel for schedule(guided)
			#pragma omp parallel for schedule(static)
			for (UINT64 i=0; i<mpi_dimension_; ++i) {
				work_[i] *= 1.0/static_cast<double>(P_[k][i]);
			}
			
			// All off-diagonal terms - no risk of data race at this point
			//#pragma omp parallel for schedule(guided)
			#pragma omp parallel for schedule(static)
			for (UINT64 i=0; i<mpi_nb_offdiag_[k]; ++i) {
				work_[mpi_local_index_friend_[k][i]] += recv_coeffs[i];
			}
		}
		
		//===========================
		// UPDATE OF LANCZOS VECTOR
		//===========================
		const double J = bond.couplingValue;
		#pragma omp parallel for schedule(static)
		for (UINT64 i=0; i<mpi_dimension_; ++i) {
			u[i] += J * work_[i];
		}
	}
	
	time(t0, "multiply");
}



} // namespace sun
