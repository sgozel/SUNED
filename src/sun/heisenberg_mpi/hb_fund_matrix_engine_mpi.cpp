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
#include "../utils/young_factor.h"
#include "../../common/time.h"
#include "../../common/mpi_utils.hpp"
#include "../../common/mpi_comm.hpp"

#ifdef SG_STORE_COLUMNS
#include "../irrep/irrep.h"
#endif

#ifdef SG_USE_VSYT
#include "../syt_usage/vsyt_usage.h"
#else
#include "../syt_usage/bsyt_usage.h"
#endif


namespace sun {

HBFundMatrixEngineMPI::HBFundMatrixEngineMPI(nlohmann::json const& inputParam)
: HBFundEngineMPI(inputParam), mvm_counter(0)
{	
	dump_matrices_ = inputParam.value("dump_matrices", false);
	
	if (dump_matrices_ == true) {
		if (!inputParam.contains("matrix_dump_folder_path")) {
			std::cerr << "Missing matrix_dump_folder_path in input .json file for HBFundMatrixEngineMPI." << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		matrix_dump_path_ = inputParam["matrix_dump_folder_path"];
		if (matrix_dump_path_.back()!='/') {
			matrix_dump_path_ += std::string("/");
		}
	}

	dump_counts_ = inputParam.value("dump_counts", false);

	if (dump_counts_ == true) {
		if (!inputParam.contains("counts_dump_folder_path")) {
			std::cerr << "Missing counts_dump_folder_path in input .json file for HBFundMatrixEngineMPI." << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		counts_dump_folder_path_ = inputParam["counts_dump_folder_path"];
		if (counts_dump_folder_path_.back()!='/') {
			counts_dump_folder_path_ += std::string("/");
		}
	}
	
	dump_runtime_ = inputParam.value("dump_runtime", false);
	n_mvm_runtime_ = inputParam.value("n_mvm_runtime", 10);

	if (dump_runtime_ == true) {
		if (!inputParam.contains("runtime_dump_folder_path")) {
			std::cerr << "Missing runtime_dump_folder_path in input .json file for HBFundMatrixEngineMPI." << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		runtime_dump_folder_path_ = inputParam["runtime_dump_folder_path"];
		if (runtime_dump_folder_path_.back()!='/') {
			runtime_dump_folder_path_ += std::string("/");
		}
	}

	{
		int max_ax = alpha_.nrows() + alpha_.ncols() - 1;
		if ((-max_ax < std::numeric_limits<typePk>::min()) || (max_ax > std::numeric_limits<typePk>::max())) {
			std::cerr << "PROBLEM : irrep has a large maximal axial distance, update typePk to support it." << std::endl;
			std::cerr << "sizeof(typePk) = " << sizeof(typePk) << std::endl;
			std::cerr << "max axial distance = " << max_ax << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}

	Y_bounds_lower_.resize(mpi_world_size_);
	
	unsigned int n_transpo = alpha_.n() - 1;
	
	if (dump_runtime_) {
		transpo_runtime_.resize(n_transpo);
		bond_runtime_.resize(lattice_.get_nbonds());
		for (unsigned int b = 0; b < lattice_.get_nbonds(); ++b) {
			bond_runtime_[b].resize(n_mvm_runtime_);
		}
	}

	P_.resize(n_transpo);
	
	local_pairs_.resize(n_transpo);
	remote_pairs_.resize(n_transpo);
	
	mpi_offdiag_nodes_remote_only_.resize(n_transpo);
	mpi_offdiag_nodes_remote_only_acc_.resize(n_transpo);
	
	for (unsigned int k=0; k<n_transpo; ++k) {
		mpi_offdiag_nodes_remote_only_[k].resize(mpi_world_size_);
		mpi_offdiag_nodes_remote_only_acc_[k].resize(mpi_world_size_);
		std::fill(mpi_offdiag_nodes_remote_only_[k].begin(), mpi_offdiag_nodes_remote_only_[k].end(), 0);
		std::fill(mpi_offdiag_nodes_remote_only_acc_[k].begin(), mpi_offdiag_nodes_remote_only_acc_[k].end(), 0);
	}
	
	mpi_nb_offdiag_.resize(n_transpo);
	mpi_nb_offdiag_local_.resize(n_transpo);
	mpi_nb_offdiag_remote_.resize(n_transpo);
	
	std::fill(mpi_nb_offdiag_.begin(), mpi_nb_offdiag_.end(), 0);
	std::fill(mpi_nb_offdiag_local_.begin(), mpi_nb_offdiag_local_.end(), 0);
	std::fill(mpi_nb_offdiag_remote_.begin(), mpi_nb_offdiag_remote_.end(), 0);
	
	mpi_local_index_base_local_.resize(n_transpo);
	mpi_local_index_friend_local_.resize(n_transpo);
	
	mpi_local_index_base_remote_.resize(n_transpo);
	mpi_local_index_friend_remote_.resize(n_transpo);
	
	if (mpi_rank_==0) {
		std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
#ifdef SG_STORE_COLUMNS
		std::cout << "WARNING: Storing SYTs as columns" << std::endl;
		std::cout << ":::::::::::::::::::::::::::::::::::::::::" << std::endl;
#endif
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
	
	// Each MPI process generates its own SYTs
	
	const UINT64 from = mpi_rank_ * mpi_bare_dimension_;
	
#ifdef SG_STORE_COLUMNS
	const Irrep alphaT(alpha_.transpose());
#ifdef SG_USE_VSYT
	std::cerr << "Currently, VSYT is unsupported on MPI application." << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
	Y_ = get_SYT<SYTel>(alphaT, from, mpi_dimension_);
#else
	Y_ = get_SYT(alphaT, from, mpi_dimension_);
#endif
#else
#ifdef SG_USE_VSYT
	std::cerr << "Currently, VSYT is unsupported on MPI application." << std::endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
	Y_ = get_SYT<SYTel>(alpha_, from, mpi_dimension_);
#else
	Y_ = get_SYT(alpha_, from, mpi_dimension_);
#endif
#endif
	
	communicate_bounds();
	
	work_.resize(mpi_dimension_);
#ifdef SG_USE_NUMA
	#pragma omp parallel for schedule(static)
	for (UINT64 i = 0; i < mpi_dimension_; ++i) {
		work_[i] = 0.0;
	}
#endif
	
	if (mpi_bare_dimension_ > std::numeric_limits<typeIndex>::max()) {
		std::cerr << "mpi_bare_dimension_ = " << mpi_bare_dimension_ << std::endl;
		std::cerr << "std::numeric_limits<typeIndex>::max() = " << std::numeric_limits<typeIndex>::max() << std::endl;
		std::cerr << "typeIndex incapable of storing mpi_bare_dimension_. Change typeIndex to a larger type." << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	time(t0, std::string("init"));
}


void HBFundMatrixEngineMPI::communicate_bounds()
{
	// pack lower bound into a vector of int
	std::vector<int8_t> local_bounds(Ns_);
	for (size_t i = 0; i < Ns_; ++i) {
		local_bounds[i] = Y_[0].get(i);
	}

	// communicate bounds
	std::vector<int8_t> all_bounds(Ns_ * mpi_world_size_);
	MPI_Allgather(
			local_bounds.data(),     // send buffer
			Ns_,                     // send count
			mpi::mpi_type<int8_t>(), // send type
			all_bounds.data(),       // recv buffer
			Ns_,                     // recv count per process
			mpi::mpi_type<int8_t>(), // recv type
			MPI_COMM_WORLD
	);
	
	// Put back vector<int8_t> into SYT format
	for (int rank = 0; rank < mpi_world_size_; ++rank) {
		for (unsigned int i = 0; i < Ns_; ++i) {
			Y_bounds_lower_[rank].set(i, all_bounds[rank * Ns_ + i]);
		}
	}
}


void HBFundMatrixEngineMPI::precise_memory_usage() const
{
	if (mpi_rank_ == 0) {
		double factor = 1e6;
		std::string units("MB");
		if (8*Y_.size() >= 1e9) {
			factor *= 1000;
			units = std::string("GB");
		}
		
#ifdef SG_USE_VSYT
		double memAllSYTs = alpha_.n() * sizeof(SYTel) * static_cast<double>(dimension_)/factor;
		double memLocalSYTs = alpha_.n() * sizeof(SYTel) * static_cast<double>(mpi_dimension_)/factor;
#else
		double memAllSYTs = sizeof(SYT_value_t) * static_cast<double>(dimension_)/factor;
		double memLocalSYTs = sizeof(SYT_value_t) * static_cast<double>(mpi_dimension_)/factor;
#endif
		
		double memLanczos = 3 * sizeof(double) * static_cast<double>(mpi_dimension_)/factor;
		double memWorkArray = sizeof(double) * static_cast<double>(mpi_dimension_)/factor;
		
		const unsigned int n_transpo = alpha_.n() - 1;
		
		double memMatrixLookups = sizeof(typePk) * n_transpo * static_cast<double>(mpi_dimension_)/factor;
		
		for (unsigned int k = 0; k < n_transpo; ++k)
		{
			double temp = 0;
			
			// mpi_local_index_base_local_[k] & mpi_local_index_friend_local_[k]
			temp += (2 * static_cast<double>(mpi_nb_offdiag_local_[k])/factor * sizeof(typeIndex));
			
			// mpi_local_index_base_remote_[k] & mpi_local_index_friend_remote_[k]
			temp += (2 * static_cast<double>(mpi_nb_offdiag_remote_[k])/factor * sizeof(typeIndex));
			
			memMatrixLookups += temp;
		}
		
		double memBufferArrays = sizeof(double) * (static_cast<double>(max_offdiag_)/factor + static_cast<double>(max_offdiag_remote_)/factor);
		
		double memTotal = memLanczos + memWorkArray + memMatrixLookups + memBufferArrays;
		
		if ((factor == 1e6) && (memTotal >= 1000))  {
			factor *= 1000;
			units = std::string("GB");
			
			memAllSYTs /= 1000;
			memLocalSYTs /= 1000;
			memLanczos /= 1000;
			memWorkArray /= 1000;
			memMatrixLookups /= 1000;
			memBufferArrays /= 1000;
			memTotal /= 1000;
		}
		
		std::ios_base::fmtflags coutflags(std::cout.flags());
		std::cout << std::fixed;
		std::cout << std::setprecision(3);
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << ":::::::: TOTAL MEMORY REQUIREMENTS ::::::::" << std::endl;
		std::cout << "::::::::     PRECISE ESTIMATE      ::::::::" << std::endl;
		std::cout << "::::::::          PER RANK         ::::::::" << std::endl;
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << "mpi_world_size_     : " << mpi_world_size_ << std::endl;
		std::cout << "---------------------------" << std::endl;
		std::cout << "Total SYTs Memory   : " << memAllSYTs << units << std::endl;
		std::cout << "Local SYTs Memory   : " << memLocalSYTs << units << std::endl;
		std::cout << "---------------------------" << std::endl;
		std::cout << "Matrix lookups      : " << memMatrixLookups << units << std::endl;
		std::cout << "3 Lanczos vectors   : " << memLanczos << units << std::endl;
		std::cout << "work_ array         : " << memWorkArray << units << std::endl;
		std::cout << "buffers             : " << memBufferArrays << units << std::endl;
		std::cout << "--------------" << std::endl;
		std::cout << "Total (per rank)    : " << memTotal << units << std::endl;
		std::cout << ":::::::::::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout.flags(coutflags);
	}
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
		
		UINT64 local_pairs = 0;
		UINT64 remote_pairs = 0;
		
		#pragma omp parallel
		{
			UINT64 local_pairs_private = 0;
			UINT64 mpi_nb_offdiag_local_k_private = 0;

			#pragma omp for schedule(guided)
			for (UINT64 i = 0; i < mpi_dimension_; ++i)
			{	
				const int rowk = Y_[i].get(k);
				const int rowkk = Y_[i].get(k+1);
				
				if (rowk == rowkk) {
					P_[k][i] = YOUNG_FACTOR;
				} else {
					const std::pair<int, int> cy = get_column_k_k_plus_one(Y_[i], k);
					if (cy.first == cy.second) {
						P_[k][i] = -YOUNG_FACTOR;
					} else {		
						SYT yfriend = Y_[i];
						yfriend.exchange(k, k+1);
						
						// Extract the rank of the friend SYT using binary search
						int rankfriend;
						{
							auto it = std::upper_bound(Y_bounds_lower_.begin(), Y_bounds_lower_.end(), yfriend);
							// it is an iterator to the first element in Y_bounds_lower which is strictly larger than yfriend
							rankfriend = static_cast<int>(std::distance(Y_bounds_lower_.begin(), it)) - 1;
						}

						if (rankfriend == -1) {
							std::cerr << "PROBLEM Y : the rank of the friend SYT could not get extracted from Y_bounds_lower_" << std::endl;
							MPI_Abort(MPI_COMM_WORLD, 1);
						}
						
						if (rankfriend == mpi_rank_) {
							// The friend SYT belongs to this rank - we can search within the local collection of SYTs
							
							// increment the private counters of local pairs
							local_pairs_private += 1;
							mpi_nb_offdiag_local_k_private += 1;
							
							// search index of element with binary search
							auto it = std::lower_bound(Y_.begin(), Y_.end(), yfriend);
							UINT64 index = it - Y_.begin(); // local index
							index += mpi_rank_ * mpi_bare_dimension_; // global index
							
							offdiag_indices[i] = index + 1; // add 1 to differentiate from value 0 used for diagonal elements
							
							// axial distance from k to k+1 in SYT Y_[i]
							// count +1 for each step made downwards or to the left
							// count -1 for each step made upwards or to the right
							const typePk ax = cy.first - rowk - cy.second + rowkk;
							
							P_[k][i] = -ax * YOUNG_FACTOR;
							
						} else {
							// register a request for an index extraction on a friend rank
							#pragma omp critical
							{
								// increment the counters of remote pairs
								remote_pairs += 1;
								mpi_nb_offdiag_remote_[k] += 1;
								mpi_offdiag_nodes_remote_only_[k][rankfriend] += 1;
								// register the request
								pending_i[rankfriend].emplace_back(i);
								pending_syt[rankfriend].emplace_back(yfriend.value());
							}
						}
					}
				}
			} // end for i (local Hilbert space)
			
			#pragma omp critical
			{
				local_pairs += local_pairs_private;
				mpi_nb_offdiag_local_[k] += mpi_nb_offdiag_local_k_private;
			}
		} // #pragma omp parallel

		// we have double-counted the local pairs, as we have counted both elements of each local pair
		if (local_pairs % 2 == 1) {
			std::cerr << "Problem : local_pairs = " << local_pairs << ", but it should be even." << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		local_pairs /= 2;
		
		local_pairs_[k] = local_pairs;
		remote_pairs_[k] = remote_pairs;
		
		if (mpi_nb_offdiag_local_[k] != 2*local_pairs_[k]) {
			std::cerr << "PROBLEM G-0 - local_pairs" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		if (mpi_nb_offdiag_remote_[k] != remote_pairs_[k]) {
			std::cerr << "PROBLEM G-1 - remote_pairs" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		// MPI EXCHANGE for remote friends
		
		if ((pending_i[mpi_rank_].size() > 0) || (pending_syt[mpi_rank_].size() > 0)) {
			std::cerr << "PROBLEM : pending_i[mpi_rank_].size() > 0 or pending_syt[mpi_rank_].size() > 0, but this rank should not communicate pending requests to itself" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
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
			mpi::mpi_type<int64_t>(),
			recv_counts.data(),
			1,
			mpi::mpi_type<int64_t>(),
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
		
		// Pack data to be sent into a contiguous buffer
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
		mpi::alltoallv(
			send_buf_syt,
			send_counts,
			send_displs,
			recv_buf_syt,
			recv_counts,
			recv_displs,
			MPI_COMM_WORLD
		);

		// search indices of received SYTs on this local collection Y_
		// and assemble the associated global indices for the MPI reply
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
		
		mpi::alltoallv(
			reply_buf,
			recv_counts,
			recv_displs,
			result_buf,
			send_counts,
			send_displs,
			MPI_COMM_WORLD
		);
		
		// Now <result_buf> contains the requested indices of the image 
		// (friend) SYTs which belong to friend ranks of local base SYTs
		
		// FILL IN PASS
		
		// we can now fill the missing information for the matrix lookups,
		// namely the data associated to the remote friend SYTs
		
		for (int r = 0; r < mpi_world_size_; ++r)
		{
			#pragma omp parallel for schedule(static)
			for (int64_t j = 0; j < send_counts[r]; ++j)
			{
				const UINT64 local_i = pending_i[r][j];
				const UINT64 global_index = result_buf[send_displs[r] + j];

				offdiag_indices[local_i] = global_index + 1; // +1 convention, same as local case

				const SYT& yi = Y_[local_i];
				const std::pair<int, int> cy = get_column_k_k_plus_one(yi, k);
				const int rowk  = yi.get(k);
				const int rowkk = yi.get(k+1);
				const typePk ax = cy.first - rowk - cy.second + rowkk; // axial distance from k to k+1
				P_[k][local_i] = -ax * YOUNG_FACTOR;
			}
		}
		
		// mpi_nb_offdiag_local_[k] = number of local off-diagonal element for transposition (k, k+1) 
		// 							  found by this process and living entirely on this process
		// 
		// mpi_nb_offdiag_remote_[k] = number of off-diagonal elements for transposition (k, k+1)
		// 							   found by this process, and having the friend SYT on a different process
		// 
		// mpi_offdiag_nodes_remote_only_[k][rank] = for transposition (k, k+1), number of remote images belonging to <rank>
		// 
		// offdiag_indices[i] = 0        for diagonal elements
		// 					  = index+1  for off-diagonal elements, with index being the global index of the image SYT
		// 
		
		if (mpi_offdiag_nodes_remote_only_[k][mpi_rank_] != 0) {
			std::cerr << "PROBLEM Z-0 - mpi_offdiag_nodes_remote_only_[k][mpi_rank_] != 0" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		for (int rank=0; rank<mpi_world_size_-1; ++rank) {
			mpi_offdiag_nodes_remote_only_acc_[k][rank+1] = mpi_offdiag_nodes_remote_only_acc_[k][rank] + mpi_offdiag_nodes_remote_only_[k][rank];
		}
		
		const UINT64 total_local_count = mpi_nb_offdiag_local_[k];
		const UINT64 total_remote_count = mpi_nb_offdiag_remote_[k];
		
		mpi_nb_offdiag_[k] = total_local_count + total_remote_count;
		
		// array which will keep track of the number of friend SYTs belonging to each rank (only remote ranks)
		std::vector<UINT64> friend_nodes_count(mpi_world_size_, 0);
		
		// arrays which will keep track of the local base index and of the friend local index for local pairs
		mpi_local_index_base_local_[k].resize(total_local_count);
		mpi_local_index_friend_local_[k].resize(total_local_count);
		
		// arrays which will keep track of the local base index on this process and of the friend index for remote pairs
		mpi_local_index_base_remote_[k].resize(total_remote_count);
		mpi_local_index_friend_remote_[k].resize(total_remote_count); // this will be the receive buffer for indices from MPI communication
		
		// 
		std::vector<typeIndex> mpi_local_index_friend_temp_remote(total_remote_count, 0);
		
		UINT64 cpt_local = 0; // will count the local pairs
		UINT64 cpt_remote = 0; // will count remote pairs
		
		for (UINT64 i=0; i<mpi_dimension_; ++i)
		{	
			if (offdiag_indices[i] > 0) {
				// extract index of image SYT
				const UINT64 index = offdiag_indices[i] - 1; // global index
				
				// extract local index of image SYT
				const UINT64 local_index_friend = mpi_local_index_from_global_index(index);
				
				// compute the MPI rank which owns this SYT
				const int rankfriend = mpi_rank_from_index(index);
				
				if (rankfriend == mpi_rank_) {
					// LOCAL pair
					
					// local index of the original (base) SYT belonging to this process
					mpi_local_index_base_local_[k][cpt_local] = static_cast<typeIndex>(i);
					
					// local index of the image (friend) SYT which also belongs to this process
					mpi_local_index_friend_local_[k][cpt_local] = static_cast<typeIndex>(local_index_friend);
					
					cpt_local += 1;
				
				} else {
					// REMOTE pair
					
					// local index in this process of the original (base) SYT
					mpi_local_index_base_remote_[k][mpi_offdiag_nodes_remote_only_acc_[k][rankfriend] + friend_nodes_count[rankfriend]] = static_cast<typeIndex>(i);
					
					// local index of image (friend) SYT in friend process
					mpi_local_index_friend_temp_remote[mpi_offdiag_nodes_remote_only_acc_[k][rankfriend] + friend_nodes_count[rankfriend]] = static_cast<typeIndex>(local_index_friend);
					
					friend_nodes_count[rankfriend] += 1;
					cpt_remote += 1;
				}
			}
		}
		
		if (friend_nodes_count[mpi_rank_] != 0) {
			std::cerr << "PROBLEM B - friend_nodes_count[mpi_rank_] = " << friend_nodes_count[mpi_rank_] << ", but should be 0" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		if (cpt_local != total_local_count) {
			std::cerr << "PROBLEM C-0 - cpt_local != total_local_count" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		if (cpt_remote != total_remote_count) {
			std::cerr << "PROBLEM C-1 - cpt_remote != total_remote_count" << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		UINT64 temp_test = 0;
		for (int rank=0; rank<mpi_world_size_; ++rank) {
			temp_test += friend_nodes_count[rank];
		}
		
		if (temp_test != total_remote_count) {
			std::cerr << "PROBLEM D - temp_test != total_remote_count" << std::endl;
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
				std::cerr << "PROBLEM X: friend_nodes_count mismatch between rank " << mpi_rank_ 
						  << " and rank " << r
						  << " friend_nodes_count[" << r << "] = " << friend_nodes_count[r]
						  << "; friend_nodes_count_received[" << r << "] = " << friend_nodes_count_received[r]
						  << std::endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}
		
		// Because of the structure of the interaction, recvcounts will be the
		// same as sendcounts, and rdispls will be the same as sdispls
		
		std::vector<int64_t>& remote_sendrecvcounts = mpi_offdiag_nodes_remote_only_[k];
		std::vector<int64_t>& remote_srdispls = mpi_offdiag_nodes_remote_only_acc_[k];
		
		mpi::alltoallv(
			mpi_local_index_friend_temp_remote,
			remote_sendrecvcounts,
			remote_srdispls,
			mpi_local_index_friend_remote_[k],
			remote_sendrecvcounts,
			remote_srdispls,
			MPI_COMM_WORLD
		);
		
		time(tk_start, std::string("Building ")+transpo_string);
		
		if (dump_matrices_ == true) {
			dump_matrix(k);
		}
		if (dump_counts_ == true) {
			if (k > 0) {
				dump_pairs_counts(k);
			}
		}
	} // for k
	
	free_basis();
	
	// count the max across remote elements (max taken over all k's)
	max_offdiag_remote_ = *std::max_element(mpi_nb_offdiag_remote_.begin(), mpi_nb_offdiag_remote_.end());
	
	// count the max of (local+remote) across all k's
	max_offdiag_ = *std::max_element(mpi_nb_offdiag_.begin(), mpi_nb_offdiag_.end());
	
	////////////////////////////////////////
	// ANALYZE LOCAL AND REMOTE PAIRS ACCROSS ALL RANKS
	////////////////////////////////////////
	
	const int n_transpo = alpha_.n() - 1;
	
	// On rank 0: receive from all ranks
	std::vector<UINT64> gathered_local;
	std::vector<UINT64> gathered_remote;

	if (mpi_rank_ == 0) {
		gathered_local.resize(n_transpo * mpi_world_size_);
		gathered_remote.resize(n_transpo * mpi_world_size_);
	}
	
	MPI_Gather(
		local_pairs_.data(),
		n_transpo,
		MPI_UINT64_T,
		gathered_local.data(),
		n_transpo,
		MPI_UINT64_T,
        0,
        MPI_COMM_WORLD
	);

	MPI_Gather(
		remote_pairs_.data(),
		n_transpo,
		MPI_UINT64_T,
		gathered_remote.data(),
		n_transpo,
		MPI_UINT64_T,
		0,
		MPI_COMM_WORLD
	);
	
	if (mpi_rank_ == 0) {
		
		// unused for the moment
		//std::vector<std::vector<UINT64>> all_local_pairs(mpi_world_size_, std::vector<UINT64>(n_transpo));
		//std::vector<std::vector<UINT64>> all_remote_pairs(mpi_world_size_, std::vector<UINT64>(n_transpo));
		
		// totals accross ranks for each transposition (k, k+1)
		std::vector<UINT64> total_local(n_transpo, 0);
		std::vector<UINT64> total_remote(n_transpo, 0);
		
		// create the matrices all_local_pairs[r][k] and all_remote_pairs[r][k] by unflattening the received arrays
		for (int r = 0; r < mpi_world_size_; ++r) {
			for (int k = 0; k < n_transpo; ++k) {
				//all_local_pairs[r][k] = gathered_local[r * n_transpo + k];
				//all_remote_pairs[r][k] = gathered_remote[r * n_transpo + k];
				
				total_local[k]  += gathered_local [r * n_transpo + k];
				total_remote[k] += gathered_remote[r * n_transpo + k];
			}
		}
		
		// correct for double-counting remote pairs due to counting each element of the pair on different ranks
		for (int k = 0; k < n_transpo; ++k) {
			if (total_remote[k] % 2 == 1) {
				std::cerr << "PROBLEM : total_remote[" << k << "] = " << total_remote[k] << " is not a multiple of 2" << std::endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			total_remote[k] /= 2;
		}
		
		// total number of pairs for each transposition (k, k+1)
		std::vector<UINT64> total_pairs(n_transpo, 0);
		for (int k = 0; k < n_transpo; ++k) {
			total_pairs[k] = total_local[k] + total_remote[k];
		}
		
		// ANALYSIS PER TRANSPOSITION
		std::cout << "==================================" << std::endl;
		std::cout << "=== Local/Remote pair analysis ===" << std::endl;
		std::cout << "==================================" << std::endl;
		for (int k = 1; k < n_transpo; ++k)
		{
			std::cout << "----------------------------------" << std::endl;
			std::cout << "Transposition (" << k << ", " << k+1 << ")" << std::endl;
			std::cout << "----------------------------------" << std::endl;
			std::cout << "dimension_         : " << dimension_ << std::endl;
			std::cout << "Total off-diags    : " << 2*total_pairs[k] << std::endl;
			std::cout << "off-diags fraction : " << std::fixed << std::setprecision(2) 
												 << 100.0 * (double)(2*total_pairs[k])/((double)(dimension_))
												 << "%" << std::endl;
			std::cout << "-----------------" << std::endl;
			std::cout << "Local pairs        : " << total_local[k] << std::endl;
			std::cout << "Remote pairs       : " << total_remote[k] << std::endl;
			std::cout << "Local fraction     : " << std::fixed << std::setprecision(2) 
											     << 100.0 * (double)(total_local[k])/((double)(total_pairs[k])) 
											     << "%" << std::endl;
		}
		std::cout << "==================================" << std::endl;
		
		// GLOBAL ANALYSIS
		UINT64 grand_total_local  = 0;
		UINT64 grand_total_remote = 0;
		UINT64 grand_total_pairs  = 0;

		for (int k = 0; k < n_transpo; ++k) {
			grand_total_local  += total_local[k];
			grand_total_remote += total_remote[k];
			grand_total_pairs  += total_pairs[k];
		}

		std::cout << "==================================" << std::endl;
		std::cout << "======== GLOBAL SUMMARY ==========" << std::endl;
		std::cout << "==================================" << std::endl;
		std::cout << "Total pairs     : " << grand_total_pairs           << std::endl;
		std::cout << "Local pairs     : " << grand_total_local           << std::endl;
		std::cout << "Remote pairs    : " << grand_total_remote          << std::endl;
		std::cout << "Local fraction  : " << std::fixed << std::setprecision(2)
										  << 100.0 * (double)grand_total_local / (double)grand_total_pairs
										  << "%" << std::endl;
		std::cout << "==================================" << std::endl;
	}
	
	////////////////////////////////////////
	
	precise_memory_usage();
	
	////////////////////////////////////////
	
	time(t_start, "build_matrix_lookups");
}


// dump local/remote counts across all ranks for transposition (k, k+1)
void HBFundMatrixEngineMPI::dump_pairs_counts(const unsigned int k)
{
	// Communicate all counts of local pairs
	std::vector<uint64_t> recv_buffer_local(mpi_world_size_, 0);

	MPI_Gather(
		&local_pairs_[k],          // send buffer
		1,                         // send count
		mpi::mpi_type<uint64_t>(), // send type
		recv_buffer_local.data(),  // receive buffer
		1,                         // receive count (per process)
		mpi::mpi_type<uint64_t>(), // send type
		0,
		MPI_COMM_WORLD
	);

	// Communicate all counts of remote pairs

	// The send buffer is:
	// mpi_offdiag_nodes_remote_only_[k] --> vector of INT64 of length mpi_world_size_
	// which contains the number of remote pairs between this rank and all other ranks

	std::vector<int64_t> recv_buffer_remote(mpi_world_size_ * mpi_world_size_, 0);

	MPI_Gather(
		mpi_offdiag_nodes_remote_only_[k].data(), // send buffer
		mpi_world_size_,           // send count
		mpi::mpi_type<int64_t>(),  // send datatype
		recv_buffer_remote.data(), // receive buffer
		mpi_world_size_,           // receive count (per process)
		mpi::mpi_type<int64_t>(),  // receive datatype
		0,                         // root (receive process)
		MPI_COMM_WORLD
	);

	if (mpi_rank_ == 0) {
		// Output to text file
		std::string filename_counts("counts_k" + std::to_string(k) + ".log");
		filename_counts = counts_dump_folder_path_ + filename_counts;
		
		std::ofstream out_counts(filename_counts, std::ios::app);
		if (!out_counts) {
			std::cerr << "Cannot open file: " + filename_counts << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		
		for (int rank1 = 0; rank1 < mpi_world_size_; ++rank1) {
			for (int rank2 = 0; rank2 < mpi_world_size_; ++rank2) {
				if (rank2  == rank1) {
					out_counts << std::left << std::setw(11) << recv_buffer_local[rank1] << "  ";
				} else {
					out_counts << std::left << std::setw(11) << recv_buffer_remote[rank1 * mpi_world_size_ + rank2] << "  ";
				}
			}
			out_counts << std::endl;
		}
	}
}


void HBFundMatrixEngineMPI::dump_runtimes() const {
	if (mpi_rank_ == 0) {
		// dump bond runtimes
		std::string filename_bond_runtimes("runtime_bonds.log");
		filename_bond_runtimes = runtime_dump_folder_path_+ filename_bond_runtimes;
		std::ofstream out_bond_runtimes(filename_bond_runtimes, std::ios::app);
		if (!out_bond_runtimes) {
			std::cerr << "Cannot open file: " + filename_bond_runtimes << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for (unsigned int b = 0; b < lattice_.get_nbonds(); ++b) {
			for (size_t i = 0; i < bond_runtime_[b].size(); ++i) {
				out_bond_runtimes << std::left << std::setw(8) 
								  << std::fixed << std::setprecision(5)
								  << bond_runtime_[b][i] << "   ";
			}
			out_bond_runtimes << std::endl;
		}

		std::vector<std::vector<double>>().swap(bond_runtime_);

		// dump transposition runtimes
		std::string filename_transpo_runtimes("runtime_transpos.log");
		filename_transpo_runtimes = runtime_dump_folder_path_+ filename_transpo_runtimes;
		std::ofstream out_transpo_runtimes(filename_transpo_runtimes, std::ios::app);
		if (!out_transpo_runtimes) {
			std::cerr << "Cannot open file: " + filename_transpo_runtimes << std::endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}

		for (unsigned int k = 0; k < alpha_.n()-1; ++k) {
			for (size_t i = 0; i < transpo_runtime_[k].size(); ++i) {
				out_transpo_runtimes << std::left << std::setw(8) 
									 << std::fixed << std::setprecision(5)
									 << transpo_runtime_[k][i] << "   ";
			}
			out_transpo_runtimes << std::endl;
		}

		std::vector<std::vector<double>>().swap(transpo_runtime_);
	}
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

	// Write mpi_nb_offdiag_local_[k]
    out.write(reinterpret_cast<const char*>(&mpi_nb_offdiag_local_[k]), sizeof(UINT64));
    
    // Write mpi_local_index_base_local_[k]
    out.write(reinterpret_cast<const char*>(mpi_local_index_base_local_[k].data()), mpi_nb_offdiag_local_[k] * sizeof(typeIndex));

    // Write mpi_local_index_friend_local_[k]
    out.write(reinterpret_cast<const char*>(mpi_local_index_friend_local_[k].data()), mpi_nb_offdiag_local_[k] * sizeof(typeIndex));
    
    // Write mpi_nb_offdiag_remote_[k]
    out.write(reinterpret_cast<const char*>(&mpi_nb_offdiag_remote_[k]), sizeof(UINT64));

	// Write mpi_local_index_base_remote_[k]
    out.write(reinterpret_cast<const char*>(mpi_local_index_base_remote_[k].data()), mpi_nb_offdiag_remote_[k] * sizeof(typeIndex));

    // Write mpi_local_index_friend_remote_[k]
    out.write(reinterpret_cast<const char*>(mpi_local_index_friend_remote_[k].data()), mpi_nb_offdiag_remote_[k] * sizeof(typeIndex));

    // Write mpi_offdiag_nodes_remote_only_[k]
    out.write(reinterpret_cast<const char*>(mpi_offdiag_nodes_remote_only_[k].data()), mpi_world_size_ * sizeof(INT64));

    // Write mpi_offdiag_nodes_remote_only_acc_[k]
    out.write(reinterpret_cast<const char*>(mpi_offdiag_nodes_remote_only_acc_[k].data()), mpi_world_size_ * sizeof(INT64));
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
	
	// Read mpi_nb_offdiag_local_[k]
    in.read(reinterpret_cast<char*>(&mpi_nb_offdiag_local_[k]), sizeof(UINT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_nb_offdiag_local_[k] from file: " << filename << std::endl;
        return false;
    }
    
    // Read mpi_local_index_base_local_[k]
    mpi_local_index_base_local_[k].resize(mpi_nb_offdiag_local_[k]);
    in.read(reinterpret_cast<char*>(mpi_local_index_base_local_[k].data()), mpi_nb_offdiag_local_[k] * sizeof(typeIndex));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_local_index_base_local_[k] from file: " << filename << std::endl;
        return false;
    }
    
    // Read mpi_local_index_friend_local_[k]
    mpi_local_index_friend_local_[k].resize(mpi_nb_offdiag_local_[k]);
    in.read(reinterpret_cast<char*>(mpi_local_index_friend_local_[k].data()), mpi_nb_offdiag_local_[k] * sizeof(typeIndex));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_local_index_friend_local_[k] from file: " << filename << std::endl;
        return false;
    }
	
    // Read mpi_nb_offdiag_remote_[k]
    in.read(reinterpret_cast<char*>(&mpi_nb_offdiag_remote_[k]), sizeof(UINT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_nb_offdiag_remote_[k] from file: " << filename << std::endl;
        return false;
    }
	
	// Read mpi_local_index_base_remote_[k]
    mpi_local_index_base_remote_[k].resize(mpi_nb_offdiag_remote_[k]);
    in.read(reinterpret_cast<char*>(mpi_local_index_base_remote_[k].data()), mpi_nb_offdiag_remote_[k] * sizeof(typeIndex));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_local_index_base_remote_[k] from file: " << filename << std::endl;
        return false;
    }

    // Read mpi_local_index_friend_remote_[k]
    mpi_local_index_friend_remote_[k].resize(mpi_nb_offdiag_remote_[k]);
    in.read(reinterpret_cast<char*>(mpi_local_index_friend_remote_[k].data()), mpi_nb_offdiag_remote_[k] * sizeof(typeIndex));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_local_index_friend_remote_[k] from file: " << filename << std::endl;
        return false;
    }
	
    // Read mpi_offdiag_nodes_remote_only_[k]
    mpi_offdiag_nodes_remote_only_[k].resize(mpi_world_size_);
    in.read(reinterpret_cast<char*>(mpi_offdiag_nodes_remote_only_[k].data()), mpi_world_size_ * sizeof(INT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_offdiag_nodes_remote_only_[k] from file: " << filename << std::endl;
        return false;
    }

    // Read mpi_offdiag_nodes_remote_only_acc_[k]
    mpi_offdiag_nodes_remote_only_acc_[k].resize(mpi_world_size_);
    in.read(reinterpret_cast<char*>(mpi_offdiag_nodes_remote_only_acc_[k].data()), mpi_world_size_ * sizeof(INT64));
    if (!in) {
        std::cerr << "Warning: failed to read mpi_offdiag_nodes_remote_only_acc_[k] from file: " << filename << std::endl;
        return false;
    }
    
    // Verify we are at end of file - no unexpected trailing data
	in.peek();
	if (!in.eof()) {
		std::cerr << "Warning: unexpected trailing data in file: " << filename << std::endl;
		return false;
	}
	
	// Deduce other attributes
	mpi_nb_offdiag_[k] = mpi_nb_offdiag_local_[k] + mpi_nb_offdiag_remote_[k];
	local_pairs_[k] = mpi_nb_offdiag_local_[k]/2;
	remote_pairs_[k] = mpi_nb_offdiag_remote_[k];
	
    return true;
}


template <class T>
std::vector<double> HBFundMatrixEngineMPI::correlations(const T& v, const unsigned int refsite) const
{
	auto t0 = std::chrono::high_resolution_clock::now();
	
	if (mpi_rank_ == 0) {
		std::cout << "::::::::::::::::::::::::::::::::::" << std::endl;
		std::cout << "Start computing correlations - reference site = " << refsite << std::endl;
		std::cout << "::::::::::::::::::::::::::::::::::" << std::endl;
	}
	
	std::vector<double> C(alpha_.n());
	C[refsite] = 1.0;
	
	buffer_.resize(max_offdiag_);
	recv_coeffs_.resize(max_offdiag_remote_);
	
	for (unsigned int t = 0; t < alpha_.n(); ++t)
	{	
		if (t != refsite) {

			// Compute C[t] = <v|(refsite, t)|v>
			
			const Transposition transpo(refsite, t);
			const std::vector<AdjacentTransposition> ops = transpo.toAdjacentTransposition();
			
			// copy eigenvector v to work_ array
			work_ = v;
			
			// apply all adjacent transpositions on work_
			apply_transpositions(ops);
			
			// Compute overlap C[t] = <v|work_>
			double overlap_loc = 0.0;
			#pragma omp parallel for reduction(+:overlap_loc) schedule(static)
			for (UINT64 i = 0; i < mpi_dimension_; ++i) {
				overlap_loc += (double)v[i] * (double)work_[i];
			}
			MPI_Allreduce(&overlap_loc, &C[t], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		}
		
		// print result to file
		if (mpi_rank_ == 0) {
			std::ios_base::fmtflags coutflags(std::cout.flags());
			std::cout << std::fixed;
			std::cout << std::setprecision(12);
			std::cout << "C[" << t << "] = " 
					  << std::setw(14) << C[t] << std::endl;
			std::cout.flags(coutflags);
		}
	} // C[t]
	
	time(t0, "correlations");
	
	return C;
}

template std::vector<double> HBFundMatrixEngineMPI::correlations<sg_vec<double>>(const sg_vec<double>&, const unsigned int) const;


void HBFundMatrixEngineMPI::apply_transpositions(const std::vector<AdjacentTransposition>& ops) const
{
	for (unsigned int j = 0; j < ops.size(); ++j)
	{	
		auto t0 = std::chrono::high_resolution_clock::now();
		
		const unsigned int k = ops[j].getk();
		
		const std::vector<int64_t>& sendrecvcounts = mpi_offdiag_nodes_remote_only_[k];
		const std::vector<int64_t>& srdispls = mpi_offdiag_nodes_remote_only_acc_[k];
		
		double* work_local_snapshot = buffer_.data();
		double* send_coeffs = buffer_.data() + mpi_nb_offdiag_local_[k];
		
		#pragma omp parallel
		{	
			// gather all remote elements to be sent to other processes
#ifdef SG_USE_NUMA
			#pragma omp for schedule(static)
#else
			#pragma omp for schedule(guided)
#endif
			for(UINT64 i = 0; i < mpi_nb_offdiag_remote_[k]; ++i) {
				const UINT64 index_base_i = mpi_local_index_base_remote_[k][i];
				const double rho = 1.0/static_cast<double>(P_[k][index_base_i]);
				send_coeffs[i] = work_[index_base_i] * std::sqrt(1.0 - rho * rho);
			}
			
			#pragma omp single nowait
			{
				//======================================
				// MPI communication of coefficients
				//======================================
				/*
				mpi::alltoallv(
					send_coeffs,
					sendrecvcounts,
					srdispls,
					recv_coeffs_.data(),
					sendrecvcounts,
					srdispls,
					MPI_COMM_WORLD
				);
				*/

				mpi::isend_irecv(
					send_coeffs,
					sendrecvcounts,
					srdispls,
					recv_coeffs_.data(),
					sendrecvcounts,
					srdispls,
					MPI_COMM_WORLD
				);
			}
		
			// gather all elements from local pairs
#ifdef SG_USE_NUMA
			#pragma omp for schedule(static)
#else
			#pragma omp for schedule(guided)
#endif
			for (UINT64 i = 0; i < mpi_nb_offdiag_local_[k]; ++i) {
				const UINT64 index_base_i = mpi_local_index_base_local_[k][i];
				const double rho = 1.0/static_cast<double>(P_[k][index_base_i]);
				work_local_snapshot[i] = work_[index_base_i] * std::sqrt(1.0 - rho * rho);
			}
		
			#pragma omp barrier
			
			//======================================
			// Perform update of <work> array
			//======================================
			
			// All diagonal terms
#ifdef SG_USE_NUMA
			#pragma omp for schedule(static)
#else
			#pragma omp for schedule(guided)
#endif
			for (UINT64 i=0; i<mpi_dimension_; ++i) {
				work_[i] *= 1.0/static_cast<double>(P_[k][i]);
			}
			
			// All local off-diagonal terms
#ifdef SG_USE_NUMA
			#pragma omp for schedule(static)
#else
			#pragma omp for schedule(guided)
#endif
			for (UINT64 i=0; i<mpi_nb_offdiag_local_[k]; ++i) {
				work_[mpi_local_index_friend_local_[k][i]] += work_local_snapshot[i];
			}
			
			// All remote off-diagonal terms
#ifdef SG_USE_NUMA
			#pragma omp for schedule(static)
#else
			#pragma omp for schedule(guided)
#endif
			for (UINT64 i=0; i<mpi_nb_offdiag_remote_[k]; ++i) {
				work_[mpi_local_index_friend_remote_[k][i]] += recv_coeffs_[i];
			}
		
		} // omp parallel section

		if ((dump_runtime_) && (mvm_counter < n_mvm_runtime_)) {
			transpo_runtime_[k].push_back(time(t0));
		}
		
	} // for j (operations in a bond)
}


void HBFundMatrixEngineMPI::multiply(const sg_vec<double> & w, sg_vec<double> & u, const double & a, const std::string & method) const
{
    if (method=="multiply_mpi_matrix_v1") {
        multiply_mpi_matrix_v1(w, u, a);
	} else {
        std::cerr << "Multiply method undefined" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
	mvm_counter += 1;

	if ((dump_runtime_) && (mvm_counter == n_mvm_runtime_)) {
		dump_runtimes();
	}
}


template <class coeff_t>
void HBFundMatrixEngineMPI::multiply_mpi_matrix_v1(const sg_vec<coeff_t>& w, sg_vec<coeff_t>& u, const double a) const
{
	// u <--- H*w - a*u
	
	auto t0 = std::chrono::high_resolution_clock::now();
	
	// u <--- -a*u
	std::for_each(u.begin(), u.end(), [a](coeff_t& el) {el*=(-a);});
	
	// pack work_local_snapshot and send_coeffs in the same buffer
	buffer_.resize(max_offdiag_);
	
	// separate buffer for received coefficients
	recv_coeffs_.resize(max_offdiag_remote_);
	
	//for (const auto& bond : lattice_.bonds)
	for (unsigned int b = 0; b < lattice_.get_nbonds(); ++b)
	{	
		auto tbond0 = std::chrono::high_resolution_clock::now();
		
		const auto& bond = lattice_.bonds[b];
		const double J = bond.couplingValue;

		if (J == 0) {
			continue;
		}

		work_ = w;
		
		apply_transpositions(bond.ops);
		
		#pragma omp parallel for schedule(static)
		for (UINT64 i = 0; i < mpi_dimension_; ++i) {
			u[i] += J * work_[i];
		}

		if ((dump_runtime_) && (mvm_counter < n_mvm_runtime_)) {
			bond_runtime_[b][mvm_counter] = time(tbond0);
		}
	}
	
	time(t0, "multiply");
}


} // namespace sun
