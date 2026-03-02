// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_IRREP_H
#define SUN_IRREP_H

#include <iostream>
#include <initializer_list>
#include <vector>
#include <algorithm>
#include <stdexcept>

namespace sun {

class Irrep
{
public:

	Irrep() { alpha_.resize(0); n_ = 0; }

	template<typename T>
    Irrep(const std::vector<T> &alpha) {
		if ((alpha.size()==0) || (alpha[0]<=0)) {
			throw std::invalid_argument("Irrep() : irrep must have at least one box.");
		}
		alpha_.resize(1, 0);
		alpha_[0] = static_cast<unsigned int>(alpha[0]);
		n_ = alpha_[0];
		for (size_t i=1; i<alpha.size(); ++i) {
			if (alpha[i]>alpha[i-1]) {
				throw std::invalid_argument("Irrep() : row lengths must be non-ascending.");
			}
			if (alpha[i]<0) {
				throw std::invalid_argument("Irrep() : row lengths must be non-negative.");
			}
			if (alpha[i]>0) {
				alpha_.push_back(static_cast<unsigned int>(alpha[i]));
				n_ += alpha_[i];
			}
		}
	}
	
	template<typename T>
    Irrep(const std::initializer_list<T> alpha) {
		for (const auto & el : alpha) {
			if (el<0) {
				throw std::invalid_argument("Irrep() : row lengths must be non-negative.");
			}
			alpha_.push_back(static_cast<unsigned int>(el));
		}
		if ((alpha_.size()==0) || (alpha_[0]==0)) {
			throw std::invalid_argument("Irrep() : irrep must have at least one box.");
		}
		n_ = alpha_[0];
		for (size_t i=1; i<alpha_.size(); ++i) {
			if (alpha_[i]>alpha_[i-1]) {
				throw std::invalid_argument("Irrep() : row lengths must be non-ascending.");
			}
			n_ += alpha_[i];
		}
		auto it = std::lower_bound(alpha_.begin(), alpha_.end(), 1, std::greater<int>());
		alpha_.erase(it, alpha_.end());
	}
	
	const unsigned int& operator[](size_t i) const {
        return alpha_[i];
    }
    
    std::vector<unsigned int> transpose() const {
		std::vector<unsigned int> alphaT(alpha_[0], 0);
		for (size_t i=0; i<alpha_.size(); ++i) {
			for (unsigned int j=0; j<alpha_[i]; ++j) {
				alphaT[j] += 1;
			}
		}
		return alphaT;
	}
    
    unsigned int n() const {return n_;}
    unsigned int nrows() const {return alpha_.size();}
    unsigned int ncols() const {
		if (n_==0) {
			return 0;
		}
		return alpha_[0];
	}
    unsigned int size() const {return alpha_.size();}
	std::vector<unsigned int> get_vector() const {return alpha_;};
	
	void print(std::ostream& os=std::cout) const {
		os << "[" << std::flush;
		size_t i=0;
		for (i=0; i<alpha_.size()-1; ++i) {
			os << alpha_[i] << ", " << std::flush;
		}
		os << alpha_[i] << "]" << std::flush;
	}
	
private:
    std::vector<unsigned int> alpha_;
    unsigned int n_;
};

} // namespace sun

#endif
