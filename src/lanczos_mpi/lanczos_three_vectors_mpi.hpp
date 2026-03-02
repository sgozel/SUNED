// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <memory>

template <class coeff_t, class type_mult>
void lanczos_step(std::vector<coeff_t> & u, 
				  std::vector<coeff_t> & v, 
				  std::vector<coeff_t> & w,
                  double & alpha, 
                  double & beta, 
                  type_mult multiply)
{
    multiply(v, w); // w <--- H*v
    
    double alpha_loc = std::inner_product(v.begin(), v.end(), w.begin(), 0.0); // v^T @ w
    MPI_Allreduce(&alpha_loc, &alpha, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    axpy<coeff_t>(w, -alpha, v); // w <--- w - alpha*v
    axpy<coeff_t>(w, -beta, u); // w <--- w - beta*u
    u = v;
    v = w;
    
    double beta_loc = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    MPI_Allreduce(&beta_loc, &beta, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    beta = std::sqrt(beta);
    
    std::for_each(v.begin(), v.end(), [beta](coeff_t& val) { val /= beta; });
}


// Lanczos diagonalization with 3 Lanczos vectors and checkpointing
template <class coeff_t, class type_mult, class type_conv>
Tmatrix lanczos(const type_mult & multiply, 
				const type_conv & converge, 
				const UINT64 dimension, 
				const LanczosParams & lp)
{
	int mpi_world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	
    Tmatrix tmat;
    
    // Generate initial Lanczos vector
    std::vector<coeff_t> v(dimension);
    lanczos_init_vector<coeff_t>(v, dimension);
    
    std::vector<coeff_t> u(dimension, 0.0);
    std::vector<coeff_t> w(dimension, 0.0);
    
    unsigned int cpt = 0;
    bool isConverged = false;
    double alpha = 0.0;
    double beta = 0.0;
    int chkpt_found = 0;        // int instead of bool to have a safe MPI communication
    int global_chkpt_found = 0; // int instead of bool to have a safe MPI communication
    
    // Try to load checkpoint
    if (lp.checkpointing && checkpoint_exists(lp.checkpoint_file)) {
		chkpt_found = 1;
	}
	
	MPI_Allreduce(&chkpt_found,        // sendbuf
				  &global_chkpt_found, // recvbuf
				  1,                   // count
				  MPI_INT,             // datatype
				  MPI_LAND,            // MPI_Op
				  MPI_COMM_WORLD);
    
    if (global_chkpt_found) {
		
		if (mpi_rank == 0) {
			std::cout << "Found existing checkpoint file. Loading ..." << std::endl;
		}
        
        LoadedCheckpoint<coeff_t> checkpoint = load_checkpoint<coeff_t>(lp.checkpoint_file);
		
        cpt = checkpoint.iteration;
        v = std::move(checkpoint.v_current);
        u = std::move(checkpoint.v_previous);
        
        tmat = Tmatrix();
        for (size_t i = 0; i < checkpoint.alpha.size(); ++i) {
            tmat.push_back(checkpoint.alpha[i], checkpoint.beta[i]);
        }
        beta = tmat.get_last_beta();
        
        if (mpi_rank == 0) {
			std::cout << "Restarting from iteration " << cpt << std::endl;
		}
    }
    
    unsigned int min_iter = (lp.min_iter < dimension - 1 ? lp.min_iter : dimension - 1);
    
    while ((cpt < min_iter) || ((isConverged == false) && (cpt < lp.max_iter)))
    {
        lanczos_step(u, v, w, alpha, beta, multiply);
        tmat.push_back(alpha, beta);
        
        if (mpi_rank == 0) {
			if (lp.logging && (cpt % lp.log_freq == 0)) {
				tmat.log(lp.logging_folder);
				tmat.log_eigval(lp.logging_folder);
			}
		}
        
        cpt += 1;
        
        if (mpi_rank == 0) {
			if ( (cpt >= min_iter) && ((lp.logging && (cpt % lp.log_freq == 0)) || (cpt % lp.conv_check_freq == 0)) ) {
				isConverged = converge(tmat);
			}
		}
		
		MPI_Bcast(
			&isConverged,
			1,
			MPI_C_BOOL,
			0, // root
			MPI_COMM_WORLD);
        
        if (lp.checkpointing && (cpt % lp.checkpoint_frequency == 0)) {
            std::vector<coeff_t> empty_vec;
            save_checkpoint(lp.checkpoint_file, 
						    cpt, 
						    tmat.get_alpha(), 
						    tmat.get_beta(), 
						    v, 
						    u, 
						    empty_vec, 
						    false, 
						    0);
        }
    }
    
    if (mpi_rank==0) {
		verify_convergence(tmat, cpt, lp, isConverged);
		
		if (lp.logging) {
			tmat.log(lp.logging_folder);
			tmat.log_eigval(lp.logging_folder);
		}
	}
    
    if ((lp.checkpointing) && (isConverged)) {
        std::remove(lp.checkpoint_file.c_str());
        std::remove((lp.checkpoint_file + ".old").c_str());
        std::remove((lp.checkpoint_file + ".new").c_str());
        
        if (mpi_rank == 0) {
			std::cout << "Converged. Checkpoint files removed." << std::endl;
		}
    }
    
    return tmat;
}


// Compute eigenvalue and first eigenvector
// On output, v_init contains the eigenvector
template <class coeff_t, class type_mult, class type_conv>
Tmatrix lanczos_eigvec(const type_mult & multiply, 
					   const type_conv & converge, 
					   const UINT64 dimension, 
					   std::vector<coeff_t>& eigvec, 
					   const LanczosParams & lp)
{
	int mpi_world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	
	eigvec.resize(0);
	eigvec.shrink_to_fit();
	
	Tmatrix tmat;
	
	// Generate initial Lanczos vector
	std::vector<coeff_t> v(dimension);
	lanczos_init_vector<coeff_t>(v, dimension);
	
    std::vector<coeff_t> u(dimension, 0.0);
    std::vector<coeff_t> w(dimension, 0.0);

    unsigned int cpt = 0;
    bool isConverged = false;

    double alpha = 0.0;
    double beta = 0.0;
    
    int chkpt_found = 0;        // int instead of bool to have a safe MPI communication
    int global_chkpt_found = 0; // int instead of bool to have a safe MPI communication
    
    bool first_pass_complete = false;
    unsigned int second_pass_cpt = 0;
    
    // Try to load checkpoint
    if (lp.checkpointing && checkpoint_exists(lp.checkpoint_file)) {
		chkpt_found = 1;
	}
	
	MPI_Allreduce(&chkpt_found,        // sendbuf
				  &global_chkpt_found, // recvbuf
				  1,                   // count
				  MPI_INT,             // datatype
				  MPI_LAND,            // MPI_Op
				  MPI_COMM_WORLD);
    
    if (global_chkpt_found) {
		
		if (mpi_rank == 0) {
			std::cout << "Found existing checkpoint file. Loading ..." << std::endl;
		}
        
        LoadedCheckpoint<coeff_t> checkpoint = load_checkpoint<coeff_t>(lp.checkpoint_file);
        
        first_pass_complete = checkpoint.first_pass_complete;
        
        if (!first_pass_complete) {
            
            if (mpi_rank == 0) {
				std::cout << "Restarting first pass (eigenvalue computation) ..." << std::endl;
			}
            cpt = checkpoint.iteration;
            v = std::move(checkpoint.v_current);
            u = std::move(checkpoint.v_previous);
            
            tmat = Tmatrix();
            for (size_t i = 0; i < checkpoint.alpha.size(); ++i) {
                tmat.push_back(checkpoint.alpha[i], checkpoint.beta[i]);
            }
            beta = tmat.get_last_beta();
            
            if (mpi_rank == 0) {
				std::cout << "Restarting from iteration " << cpt << std::endl;
			}
            
        } else {
            
            if (mpi_rank == 0) {
				std::cout << "First pass complete. Restarting second pass (eigenvector computation) ... " << std::endl;
			}
            
            cpt = checkpoint.iteration;
            second_pass_cpt = checkpoint.second_pass_iteration;
            
            tmat = Tmatrix();
            for (size_t i = 0; i < checkpoint.alpha.size(); ++i) {
                tmat.push_back(checkpoint.alpha[i], checkpoint.beta[i]);
            }
            beta = tmat.get_beta()[second_pass_cpt-2];
            
            // Restore second pass state
            v = std::move(checkpoint.v_current);
            u = std::move(checkpoint.v_previous);
            
            eigvec.resize(dimension);
            eigvec = std::move(checkpoint.eigvec);
            
            if (mpi_rank == 0) {
				std::cout << "Restarting second pass from iteration " << second_pass_cpt << std::endl;
			}
        }
    }
    
    // ========== FIRST PASS: Compute eigenvalue ==========
    if (!first_pass_complete) {
		
        unsigned int min_iter = (lp.min_iter < dimension - 1 ? lp.min_iter : dimension - 1);
        
        while ((cpt < min_iter) || ((isConverged == false) && (cpt < lp.max_iter)))
        {
            lanczos_step(u, v, w, alpha, beta, multiply);
            tmat.push_back(alpha, beta);
            
            if (mpi_rank==0) {
				if (lp.logging & (cpt % lp.log_freq == 0)) {
					tmat.log(lp.logging_folder);
					tmat.log_eigval(lp.logging_folder);
				}
			}
            
            cpt += 1;
            
            if (mpi_rank == 0) {
				if ( (cpt >= min_iter) && ((lp.logging && (cpt % lp.log_freq == 0)) || (cpt % lp.conv_check_freq == 0)) ) {
					isConverged = converge(tmat);
				}
			}
            
            MPI_Bcast(
				&isConverged, 	// send data
				1, 				// count
				MPI_C_BOOL,  	// MPI_Datatype
				0, 				// root
				MPI_COMM_WORLD);
            
            if (lp.checkpointing && (cpt % lp.checkpoint_frequency == 0)) {
                const std::vector<coeff_t> emtpy_vec(0);
                save_checkpoint(lp.checkpoint_file, 
 							    cpt, 
							    tmat.get_alpha(), 
							    tmat.get_beta(), 
							    v, 
							    u, 
							    emtpy_vec, 
							    false, 
							    0);
            }
        }
        
        if (mpi_rank==0) {
			verify_convergence(tmat, cpt, lp, isConverged);
			
			if (lp.logging) {
				tmat.log(lp.logging_folder);
				tmat.log_eigval(lp.logging_folder);
			}
		}
        first_pass_complete = true;
    }
    
    // ========== SECOND PASS: Build eigenvector ==========
    std::cout << "Second pass -- build eigenvector ..." << std::endl;
    
    std::vector<std::vector<double>> eigvecs = tmat.eigenvectors();
    const auto& gs = eigvecs[0];
    
    // Initialize second pass (only if not restarting from second pass checkpoint)
    if (second_pass_cpt == 0) {
        
        lanczos_init_vector<coeff_t>(v, dimension);
        std::fill(u.begin(), u.end(), 0.0);
        std::fill(w.begin(), w.end(), 0.0);
        
        eigvec.resize(dimension);
        eigvec.shrink_to_fit();
        eigvec = v;
        
        const double gs0 = gs[0];
        std::for_each(eigvec.begin(), eigvec.end(), [gs0](coeff_t& el) { el *= gs0; });
        
        beta = 0.0; // not stricly necessary, because u has been re-initialized to 0
        second_pass_cpt += 1;
    }
    
    // Continue or start second pass
    for (unsigned int j = second_pass_cpt; j < cpt; ++j) {
        
        lanczos_step(u, v, w, alpha, beta, multiply);
        
        axpy<coeff_t>(eigvec, gs[j], v);
        
        second_pass_cpt = j + 1;
        
        if (lp.checkpointing && (second_pass_cpt % lp.checkpoint_frequency == 0)) {
            save_checkpoint(lp.checkpoint_file, 
						    cpt, 
						    tmat.get_alpha(), 
						    tmat.get_beta(), 
						    v, 
						    u, 
						    eigvec, 
						    true, 
						    second_pass_cpt);
        }
    }
    
    double norm_loc = std::sqrt(std::inner_product(eigvec.begin(), eigvec.end(), eigvec.begin(), 0.0));
    double norm = 0.0;
    MPI_Allreduce(&norm_loc, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    std::for_each(eigvec.begin(), eigvec.end(), [norm](coeff_t& el) { el /= norm; });
    
    if (lp.dump_eigvec == true) {
		dump_eigvec(eigvec, 0, lp);
	}
    
    if ((lp.checkpointing) && (isConverged)) {
        std::remove(lp.checkpoint_file.c_str());
        std::remove((lp.checkpoint_file + ".old").c_str());
        std::remove((lp.checkpoint_file + ".new").c_str());
        if (mpi_rank == 0) {
			std::cout << "Converged. Checkpoint files removed." << std::endl;
		}
    }
    
	return tmat;
}
