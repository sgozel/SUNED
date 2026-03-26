// Copyright 2026 Samuel GOZEL, GNU GPLv3


template<class coeff_t>
void lanczos_init_vector(sg_vec<coeff_t>& v, const UINT64 dimension, const unsigned int seed)
{	
	v.resize(dimension);
	v.shrink_to_fit();
    
    // random vector
	std::mt19937 gen(seed);
	std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::generate(v.begin(), v.end(), [&]() {return dist(gen);});
    
    // testing : uniform vector
    //std::generate(v.begin(), v.end(), [&]() {return 1.0;});
    
    double norm = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    
#ifdef SG_USE_MPI
    double global_norm;
    MPI_Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    norm = global_norm;
#endif
    norm = std::sqrt(norm);
    std::for_each(v.begin(), v.end(), [norm](coeff_t& el) { el/=norm; });
}


template<class coeff_t>
void numa_lanczos_init_vector(sg_vec<coeff_t>& v, const UINT64 dimension, const unsigned int seed)
{	
	v.resize(dimension);
	v.shrink_to_fit();
    
    // Parallel random fill — each thread gets its own seed
    // This ensures first-touch distributes pages across NUMA nodes
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        std::mt19937 gen(seed + tid);
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        
        #pragma omp for schedule(static)
        for (UINT64 i = 0; i < dimension; ++i) {
            v[i] = dist(gen);
        }
    }
    
    double norm = 0.0;
    #pragma omp parallel for reduction(+:norm) schedule(static)
    for (UINT64 i = 0; i < dimension; ++i) {
        norm += (double)v[i] * (double)v[i];
    }
    
#ifdef SG_USE_MPI
    double global_norm;
    MPI_Allreduce(&norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    norm = global_norm;
#endif
    norm = std::sqrt(norm);
    
    #pragma omp parallel for schedule(static)
    for (UINT64 i = 0; i < dimension; ++i) {
		v[i] /= norm;
	}
}


std::vector<double> residual(Tmatrix & tmat, unsigned int k)
{
	const unsigned int n = tmat.size();
	if (n == 1) {
		throw std::runtime_error("Impossible to compute residual of 1-dimensional Tmatrix");
	}
	k = std::min(k, n-1);
	Tmatrix tmat_prev = tmat;
	tmat_prev.pop_back();
	std::vector<std::vector<double>> eigvecs_prev = tmat_prev.eigenvectors();
	
	std::vector<double> residuals(k, 0.0);
	for (unsigned int i = 0; i < k; ++i) {
		residuals[i] = std::abs(tmat.get_last_beta() * eigvecs_prev[i][tmat_prev.size()-1]);
	}
	return residuals;
}


std::vector<double> ritz_value_stabilization(Tmatrix & tmat, unsigned int k)
{
	const unsigned int n = tmat.size();
	if (n == 1) {
		throw std::runtime_error("Impossible to compute ritz_value_stabilization of 1-dimensional Tmatrix");
	}
	k = std::min(k, n-1);
    std::vector<double> eigvals = tmat.eigenvalues();
    Tmatrix tmat_prev = tmat;
    tmat_prev.pop_back();
    std::vector<double> eigvals_prev = tmat_prev.eigenvalues();
    
    std::vector<double> rvs(k, 0.0);
    for (unsigned int i = 0; i < k; ++i) {
		rvs[i] = std::abs(eigvals[i] - eigvals_prev[i]) / std::abs(eigvals[i]);
	}
    return rvs;
}


bool convergence(Tmatrix & tmat, const LanczosParams & lp)
{
    bool isConverged = false;
    const unsigned int n = tmat.size();

	if (n == 1) {
		return isConverged;
	}

	const unsigned int k = std::min(lp.k, n-1);

	std::vector<double> residuals = residual(tmat, k);
	std::vector<double> rvs = ritz_value_stabilization(tmat, k);

	if (lp.logging == true) {
		// log residuals
		std::string filename_residuals = lp.logging_folder + "residuals.log";
		std::ofstream out_residuals(filename_residuals, std::ios::app);
		if (!out_residuals) {
			throw std::runtime_error("Cannot open file: " + filename_residuals);
		}
		out_residuals << n << ": " << std::flush;
		out_residuals << std::fixed << std::setprecision(12) << std::scientific;
		for (size_t i = 0; i < residuals.size(); ++i) {
			out_residuals << residuals[i] << std::flush << " ";
		}
		out_residuals << std::endl;
		// log rvs
		std::string filename_rvs = lp.logging_folder + "rvs.log";
		std::ofstream out_rvs(filename_rvs, std::ios::app);
		if (!out_rvs) {
			throw std::runtime_error("Cannot open file: " + filename_rvs);
		}
		out_rvs << n << ": " << std::flush;
		out_rvs << std::fixed << std::setprecision(12) << std::scientific;
		for (size_t i = 0; i < rvs.size(); ++i) {
			out_rvs << rvs[i] << std::flush << " ";
		}
		out_rvs << std::endl;
	}

	if (k >= lp.k) {
		double tol_residual = lp.tol_residual;
    	double tol_ritz = lp.tol_ritz;
        
        bool b1 = std::all_of(residuals.begin(), residuals.end(), [&tol_residual](const auto& el) { return el<tol_residual; });
        bool b2 = std::all_of(rvs.begin(), rvs.end(), [tol_ritz](const auto& el) {return el<tol_ritz;});
        if (b1 || b2) {
			isConverged = true;
		}
    }
    return isConverged;
}


void verify_convergence(Tmatrix & tmat, const unsigned int cpt, const LanczosParams & lp, const bool isConverged)
{
	std::cout << "Performed " << cpt << " iterations" << std::endl;

    if ( (cpt == lp.max_iter) && (isConverged == false) ) {
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "Warning: in Lanczos: reached maximum number of iterations (max_iter=" << lp.max_iter << ")." << std::endl;
        if (tmat.size() != lp.max_iter) {
			throw std::runtime_error("Problem: T matrix is not of dimension max_iter.");
		}
	}
	
	std::cout << "Tmatrix is of dimension : " << tmat.size() << std::endl;
	
	const unsigned int k = std::min(cpt, std::max(lp.k+1, static_cast<unsigned int>(10)));
	std::vector<double> residuals = residual(tmat, k);
	std::vector<double> rvs = ritz_value_stabilization(tmat, k);
	std::vector<double> eigvals = tmat.eigenvalues();
	
	std::ios_base::fmtflags coutflags(std::cout.flags());
	for (size_t t = 0; t < residuals.size(); ++t) {
		std::cout << "Eigenvalue[" << t << "] = "
				  << std::setprecision(16) << std::defaultfloat << eigvals[t] 
				  << " (r="
				  << std::scientific << std::setprecision(4) << residuals[t] 
				  << "; rvs=" 
				  << std::scientific << std::setprecision(4) << rvs[t] 
				  << ")" << std::endl;
	}
	std::cout.flags(coutflags);
	std::cout << "--------------------------------------" << std::endl;
}


template<typename E, typename coeff_t, class Alloc>
void dump_eigpair(const E energy, 
				  const std::vector<coeff_t, Alloc>& eigvec, 
				  const unsigned int index, 
				  const LanczosParams & lp)
{
#ifdef SG_USE_MPI
	int mpi_world_size;
	int mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	std::string filename = lp.eigvec_folder 
						 + std::string("eigvec_") + std::to_string(index) 
						 + "_ws" + std::to_string(mpi_world_size)
						 + "_rank" + std::to_string(mpi_rank) 
						 + "_seed" + std::to_string(lp.seed) 
						 + ".bin";
#else
	std::string filename = lp.eigvec_folder 
					     + std::string("eigvec_") + std::to_string(index) 
					     + "_seed" + std::to_string(lp.seed) 
					     + ".bin";
#endif
	
	std::ofstream out(filename, std::ios::binary);
	if (!out) {
		throw std::runtime_error("Cannot open eigvec file: " + filename);
	}
	
	// Write dimension
	const UINT64 dim = eigvec.size();
	out.write(reinterpret_cast<const char*>(&dim), sizeof(UINT64));
	
	// Write energy
	out.write(reinterpret_cast<const char*>(&energy), sizeof(E));
	
	// Write eigenvector
	out.write(reinterpret_cast<const char*>(eigvec.data()), dim * sizeof(coeff_t));
    
    out.flush();
	if (!out.good()) {
		throw std::runtime_error("Error writing eigvec to file");
	}
}


template<typename E, typename coeff_t, class Alloc>
bool load_eigpair(std::pair<E, std::vector<coeff_t, Alloc>>& eigpair, 
				  const unsigned int index, 
				  const LanczosParams & lp)
{
#ifdef SG_USE_MPI
	int mpi_world_size;
	int mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	std::string filename = lp.eigvec_folder 
						 + std::string("eigvec_") + std::to_string(index) 
						 + "_ws" + std::to_string(mpi_world_size)
						 + "_rank" + std::to_string(mpi_rank) 
						 + "_seed" + std::to_string(lp.seed) 
						 + ".bin";
#else
	std::string filename = lp.eigvec_folder 
					     + std::string("eigvec_") + std::to_string(index) 
					     + "_seed" + std::to_string(lp.seed) 
					     + ".bin";
#endif
	
	std::ifstream in(filename, std::ios::binary);
	if (!in) {
		throw std::runtime_error("Cannot open eigvec file: " + filename);
	}
	
	// Read eigvec dimension
    UINT64 dim = 0;
    in.read(reinterpret_cast<char*>(&dim), sizeof(UINT64));
    if (!in) {
        std::cerr << "Warning: failed to read dimension in file: " << filename << std::endl;
        return false;
    }
    
    // Read energy
    in.read(reinterpret_cast<char*>(&eigpair.first), sizeof(E));
    if (!in) {
        std::cerr << "Warning: failed to read energy from file: " << filename << std::endl;
        return false;
    }
    
    // Read eigenvector
    eigpair.second.resize(dim);
    in.read(reinterpret_cast<char*>(eigpair.second.data()), dim * sizeof(coeff_t));
    if (!in) {
        std::cerr << "Warning: failed to read eigenvector from file: " << filename << std::endl;
        return false;
    }
    
    // Verify we are at end of file - no unexpected trailing data
	in.peek();
	if (!in.eof()) {
		std::cerr << "Warning: unexpected trailing data in file: " << filename << std::endl;
		return false;
	}
	
	return true;
}


template<class coeff_t>
void axpy(sg_vec<coeff_t>& y, coeff_t a, const sg_vec<coeff_t>& x) {
    for (size_t i=0; i<y.size(); ++i) {
        y[i] += a * x[i];
    }
}


template<class coeff_t>
void numa_axpy(sg_vec<coeff_t>& y, coeff_t a, const sg_vec<coeff_t>& x) {
    const UINT64 dimension = y.size();
    #pragma omp parallel for schedule(static)
    for (UINT64 i = 0; i < dimension; ++i) {
        y[i] += a * x[i];
	}
}
