// Copyright 2026 Samuel GOZEL, GNU GPLv3

/* COMPARE THIS FILE WITH THE NORMAL LANCZOS_CHECKPOINTING.HPP !!!
 * 
 * BOTH SHOULD BE THE SAME
 * 
 * IT IS IS THE CASE, REPLACE THIS ONE BY THE NORMAL LANCZOS_CHECKPOINTING !!!
 * 
 */


template<class coeff_t>
struct LoadedCheckpoint {
    unsigned int iteration;
    std::vector<double> alpha;
    std::vector<double> beta;
    sg_vec<coeff_t> v_current;
    sg_vec<coeff_t> v_previous;
    bool first_pass_complete;
    unsigned int second_pass_iteration;
    sg_vec<coeff_t> eigvec;
};


inline bool checkpoint_exists(const std::string& filename) {
    std::ifstream in(filename, std::ios::binary);
    return in.good();
}


template<class coeff_t, class Alloc>
void save_checkpoint(const std::string& filename, 
					 unsigned int iteration, 
					 const std::vector<double>& alpha, 
					 const std::vector<double>& beta, 
					 const std::vector<coeff_t, Alloc>& v_current, 
					 const std::vector<coeff_t, Alloc>& v_previous, 
					 const std::vector<coeff_t, Alloc>& eigvec, 
					 bool first_pass_complete, 
					 unsigned int second_pass_iteration)
{
	std::cout << "Start saving checkpoint ..." << std::endl;
	
    std::string temp_file = filename + ".new";
    std::string backup_file = filename + ".old";
    
    // Write to temporary file
    {
        std::ofstream out(temp_file, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Cannot open temporary checkpoint file: " + temp_file);
        }
        
        // Write first pass iteration count
        out.write(reinterpret_cast<const char*>(&iteration), sizeof(iteration));
        
        // Write pass status information
        out.write(reinterpret_cast<const char*>(&first_pass_complete), sizeof(first_pass_complete));
        
        // Write second pass iteration count
        out.write(reinterpret_cast<const char*>(&second_pass_iteration), sizeof(second_pass_iteration));
        
        // Write size of alpha
        const size_t alpha_size = alpha.size();
        if (alpha_size != iteration) {
			throw std::runtime_error("alpha.size() != iteration.");
		}
        out.write(reinterpret_cast<const char*>(&alpha_size), sizeof(alpha_size));
        
        // Write alpha
        out.write(reinterpret_cast<const char*>(alpha.data()), alpha_size * sizeof(double));
        
        // Write beta
        if (beta.size() != alpha_size) {
			throw std::runtime_error("alpha and beta do not have the same size.");
		}
        //out.write(reinterpret_cast<const char*>(&beta_size), sizeof(beta_size));
        out.write(reinterpret_cast<const char*>(beta.data()), alpha_size * sizeof(double));
        
        // Write vector size (Hilbert space dimension)
        const UINT64 dimension = v_current.size();
        out.write(reinterpret_cast<const char*>(&dimension), sizeof(dimension));
        
        // Write v_current
        out.write(reinterpret_cast<const char*>(v_current.data()), dimension * sizeof(coeff_t));
        
        // Write v_previous
        if (v_previous.size() != dimension) {
			throw std::runtime_error("v_previous does not have the same dimension as v_current.");
		}
        out.write(reinterpret_cast<const char*>(v_previous.data()), dimension * sizeof(coeff_t));
        
        // Write eigvec dimension
        const UINT64 eigvec_dimension = eigvec.size();  // 0 in first pass
		out.write(reinterpret_cast<const char*>(&eigvec_dimension), sizeof(eigvec_dimension));
		
		// Write eigvec
		if (eigvec_dimension > 0) {
			if (eigvec_dimension != dimension) {
				throw std::runtime_error("non-empty eigvec does not have the same dimension as v_current.");
			}
			out.write(reinterpret_cast<const char*>(eigvec.data()), 
					  eigvec_dimension * sizeof(coeff_t));
		}
        
        out.flush();
        if (!out.good()) {
            throw std::runtime_error("Error writing checkpoint");
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Atomic rename
    if (checkpoint_exists(filename)) {
        std::remove(backup_file.c_str());
        if (std::rename(filename.c_str(), backup_file.c_str()) != 0) {
            throw std::runtime_error("Failed to rename old checkpoint");
        }
    }
    
    if (std::rename(temp_file.c_str(), filename.c_str()) != 0) {
        throw std::runtime_error("Failed to rename new checkpoint");
    }
    
    std::remove(backup_file.c_str());
    
    if (first_pass_complete) {
        std::cout << "Checkpoint saved: second pass iteration " << second_pass_iteration 
                  << " (of " << iteration << " total)" << std::endl;
    } else {
        std::cout << "Checkpoint saved: first pass iteration " << iteration << std::endl;
    }
}


template<class coeff_t>
void load_checkpoint(const std::string& filename, LoadedCheckpoint<coeff_t>& checkpoint)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Cannot open checkpoint file: " + filename);
    }
    
    // Read first pass iteration count
    in.read(reinterpret_cast<char*>(&checkpoint.iteration), sizeof(checkpoint.iteration));
    
    // Read pass status information
    in.read(reinterpret_cast<char*>(&checkpoint.first_pass_complete), sizeof(checkpoint.first_pass_complete));
    
    // Read second pass iteration count
    in.read(reinterpret_cast<char*>(&checkpoint.second_pass_iteration), sizeof(checkpoint.second_pass_iteration));
    
    // Read size of alpha
    size_t alpha_size;
    in.read(reinterpret_cast<char*>(&alpha_size), sizeof(alpha_size));
    
    // Read alpha
    checkpoint.alpha.resize(alpha_size);
    in.read(reinterpret_cast<char*>(checkpoint.alpha.data()), alpha_size * sizeof(double));
    
    // Read beta
    checkpoint.beta.resize(alpha_size);
    in.read(reinterpret_cast<char*>(checkpoint.beta.data()), alpha_size * sizeof(double));
    
    // Read vector size (Hilbert space dimension)
    UINT64 dimension;
    in.read(reinterpret_cast<char*>(&dimension), sizeof(dimension));
    
    // Read v_current
    checkpoint.v_current.resize(dimension);
    in.read(reinterpret_cast<char*>(checkpoint.v_current.data()), dimension * sizeof(coeff_t));
    // Re-touch to restore NUMA layout
	{
		sg_vec<coeff_t> tmp(dimension);
		#pragma omp parallel for schedule(static)
		for (UINT64 i = 0; i < dimension; ++i) {
			tmp[i] = checkpoint.v_current[i];
		}
		std::swap(checkpoint.v_current, tmp);
	}
    
    // Read v_previous
    checkpoint.v_previous.resize(dimension);
    in.read(reinterpret_cast<char*>(checkpoint.v_previous.data()), dimension * sizeof(coeff_t));
    // Re-touch to restore NUMA layout
	{
		sg_vec<coeff_t> tmp(dimension);
		#pragma omp parallel for schedule(static)
		for (UINT64 i = 0; i < dimension; ++i) {
			tmp[i] = checkpoint.v_previous[i];
		}
		std::swap(checkpoint.v_previous, tmp);
	}
    
    // Read eigvec dimension
    UINT64 eigvec_dimension;
	in.read(reinterpret_cast<char*>(&eigvec_dimension), sizeof(eigvec_dimension));
	checkpoint.eigvec.resize(eigvec_dimension);
	
	// Read eigvec
	if (eigvec_dimension > 0) {
		in.read(reinterpret_cast<char*>(checkpoint.eigvec.data()), 
				eigvec_dimension * sizeof(coeff_t));
		// Re-touch to restore NUMA layout
		{
			sg_vec<coeff_t> tmp(eigvec_dimension);
			#pragma omp parallel for schedule(static)
			for (UINT64 i = 0; i < eigvec_dimension; ++i) {
				tmp[i] = checkpoint.eigvec[i];
			}
			std::swap(checkpoint.eigvec, tmp);
		}
	}
    
    if (!in.good() && !in.eof()) {
        throw std::runtime_error("Error reading checkpoint file");
    }
}


template<class coeff_t>
LoadedCheckpoint<coeff_t> load_checkpoint(const std::string& filename)
{
    std::string backup_file = filename + ".old";
    LoadedCheckpoint<coeff_t> checkpoint;
    
    bool loaded = false;
    try {
        load_checkpoint(filename, checkpoint);
        std::cout << "Successfully loaded checkpoint" << std::endl;
        loaded = true;
    } catch (const std::exception& e1) {
        std::cerr << "Failed to load main checkpoint: " << e1.what() << std::endl;
        if (checkpoint_exists(backup_file)) {
            std::cout << "Attempting to load backup checkpoint ..." << std::endl;
            try {
                load_checkpoint(backup_file, checkpoint);
                loaded = true;
                std::cout << "Successfully loaded backup checkpoint" << std::endl;
                std::remove(filename.c_str());
                std::rename(backup_file.c_str(), filename.c_str());
            } catch (const std::exception& e2) {
                std::cerr << "Failed to load backup checkpoint: " << e2.what() << std::endl;
            }
        }
    }
    
    if (!loaded) {
        throw std::runtime_error("Failed to load any valid checkpoint file");
    }
    
    return checkpoint;
}
