// Copyright 2026 Samuel GOZEL, GNU GPLv3


// Perform one Lanczos step with 2 Lanczos vectors
template <class coeff_t, class type_mult>
inline void lanczos_step_two_vectors(std::vector<coeff_t> & v, 
									 std::vector<coeff_t> & w, 
									 double & alpha, 
									 double & beta, 
									 type_mult multiply)
{
	std::for_each(w.begin(), w.end(), [beta] (coeff_t& el) {el/=beta;});	
	std::for_each(v.begin(), v.end(), [beta] (coeff_t& el) {el*=(-beta);});
	std::swap(v, w);
	multiply(v, w, -1.0); // w <--- H*v + w
	alpha = std::inner_product(v.begin(), v.end(), w.begin(), 0.0);
	axpy<coeff_t>(w, -alpha, v); // w <--- w - alpha*v
	beta = std::sqrt(std::inner_product(w.begin(), w.end(), w.begin(), 0.0));
}


// Lanczos diagonalization with 2 Lanczos vectors -- no output eigenvector
template <class coeff_t, class type_mult, class type_conv>
Tmatrix lanczos_two_vectors(const type_mult & multiply, 
							const type_conv & converge, 
							const UINT64 dimension, 
							const LanczosParams & lp)
{
	Tmatrix tmat;
    
    // Generate initial Lanczos vector
    std::vector<coeff_t> w(dimension);
    lanczos_init_vector<coeff_t>(w, dimension, lp.seed);
    std::vector<coeff_t> v(dimension, 0.0);

    unsigned int cpt = 0;
    bool isConverged = false;

    double alpha = 0.0;
    double beta = 1.0;
    
    unsigned int min_iter = (lp.min_iter<dimension-1 ? lp.min_iter : dimension-1);
	
    while ((cpt<min_iter) || ((isConverged==false) && (cpt<lp.max_iter)))
    {
		lanczos_step_two_vectors(v, w, alpha, beta, multiply);
		tmat.push_back(alpha, beta);
		
		if (lp.logging && (cpt % lp.log_freq == 0)) {
			tmat.log(lp.logging_folder);
			tmat.log_eigval(lp.logging_folder);
		}
		
		cpt += 1;
		
		if ( (cpt >= min_iter) && ((lp.logging && (cpt % lp.log_freq == 0)) || (cpt % lp.conv_check_freq == 0)) ) {
			isConverged = converge(tmat);
		}
	}
	
	verify_convergence(tmat, cpt, lp, isConverged);
    if (lp.logging) {
		tmat.log(lp.logging_folder);
		tmat.log_eigval(lp.logging_folder);
	}
    
    return tmat;
}



// Lanczos diagonalization with 2 Lanczos vectors -- with output eigenvector
template <class coeff_t, class type_mult, class type_conv>
Tmatrix lanczos_two_vectors_eigvec(const type_mult & multiply, 
								   const type_conv & converge, 
								   const UINT64 dimension, 
								   std::vector<coeff_t> & eigvec, 
								   const LanczosParams & lp)
{
	Tmatrix tmat;
	
	std::vector<coeff_t> w(dimension);
    lanczos_init_vector<coeff_t>(w, dimension, lp.seed);
    std::vector<coeff_t> v(dimension, 0.0);

    unsigned int cpt = 0;
    bool isConverged = false;

    double alpha = 0.0;
    double beta = 1.0;
    
    //bool first_pass_complete = false; // to be used with checkpointing
    unsigned int second_pass_cpt = 0;
    
    unsigned int min_iter = (lp.min_iter<dimension-1 ? lp.min_iter : dimension-1);
	
	// ========== FIRST PASS: Compute eigenvalue ==========
	
    while ((cpt<min_iter) || ((isConverged==false) && (cpt<lp.max_iter)))
    {
		lanczos_step_two_vectors(v, w, alpha, beta, multiply);
		tmat.push_back(alpha, beta);
		
		if (lp.logging && (cpt % lp.log_freq == 0)) {
			tmat.log(lp.logging_folder);
			tmat.log_eigval(lp.logging_folder);
		}
		
		cpt += 1;
		
		if ( (cpt >= min_iter) && ((lp.logging && (cpt % lp.log_freq == 0)) || (cpt % lp.conv_check_freq == 0)) ) {
			isConverged = converge(tmat);
		}
	}
	
    verify_convergence(tmat, cpt, lp, isConverged);
    if (lp.logging) {
		tmat.log(lp.logging_folder);
		tmat.log_eigval(lp.logging_folder);
	}
	
	//first_pass_complete = true;
    
	// ========== SECOND PASS: Build eigenvector ==========
    std::cout << "Second pass -- build eigenvector ..." << std::endl;
    
    std::vector<std::vector<double>> eigvecs = tmat.eigenvectors();
    const auto& gs = eigvecs[0];
    
    lanczos_init_vector<coeff_t>(w, dimension, lp.seed);
    std::fill(v.begin(), v.end(), 0.0);
    eigvec.resize(dimension);
    eigvec.shrink_to_fit();
    eigvec = w;
	
	const double gs0 = gs[0];
    std::for_each(eigvec.begin(), eigvec.end(), [gs0] (coeff_t& el) { el *= gs0; });
    beta = 1.0; // strictly necessary
    second_pass_cpt += 1;
    
    for (unsigned int j=second_pass_cpt; j<cpt; ++j) {
		lanczos_step_two_vectors(v, w, alpha, beta, multiply);
		axpy<coeff_t>(eigvec, gs[j]/beta, w);
	}
	
	double norm = std::sqrt(std::inner_product(eigvec.begin(), eigvec.end(), eigvec.begin(), 0.0));
	std::for_each(eigvec.begin(), eigvec.end(), [norm] (coeff_t& el) {el/=norm;});
    
    if (lp.dump_eigvec == true) {
		dump_eigvec<coeff_t>(eigvec, 0, lp);
	}
    
    return tmat;
}
