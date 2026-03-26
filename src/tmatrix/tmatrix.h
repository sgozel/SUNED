// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef TMATRIX_H
#define TMATRIX_H

#include <vector>
#include <string>
#include <map>
#include <utility>

#ifdef SG_USE_MPI
namespace lanczosmpi {
#else
namespace lanczos {
#endif

class Tmatrix
{
public:
    Tmatrix();
    Tmatrix(const std::vector<double>&, const std::vector<double> &);

    void push_back(const double &, const double &);
    void pop_back();

    unsigned int size() const {return n_;};
	
    std::vector<double> get_alpha() const {return alpha_;};
    std::vector<double> get_beta() const {return beta_;};
    double get_last_beta() const {return beta_[n_-1];};

    std::vector<double> eigenvalues();
    std::vector<std::vector<double>> eigenvectors();
    std::pair<std::vector<double>, std::vector<std::vector<double>>> eig();
    
    void log(const std::string& folder = std::string("")) const;
    void log_eigvals(unsigned int k = 1, const std::string& folder = std::string("")) const;
    void log_full(const std::string& folder = std::string("")) const;

protected:
    std::vector<double> alpha_;
    std::vector<double> beta_;
    unsigned int n_;
    
    std::map<unsigned int, std::vector<double>> eigvals_;
};

}

#endif
