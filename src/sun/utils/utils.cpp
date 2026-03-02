// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "utils.h"

#include <numeric>
#include <stdexcept>

namespace sun {

void reduce_by_divide(std::vector<unsigned int>& numvec, std::vector<unsigned int>& denomvec) {
	
	std::sort(numvec.begin(), numvec.end());
    std::sort(denomvec.begin(), denomvec.end());
	
	// Keep only values > 1
	auto f = [] (const auto & el) -> bool {return el==1;};
	
	std::erase_if(numvec, f);
	std::erase_if(denomvec, f);
    
    if (denomvec.empty()) {
        if (numvec.empty()) {
            numvec = {1};
            denomvec = {1};
            return;
        } else {
            denomvec = {1};
            return;
        }
    } else {
        if (numvec.empty()) {
            throw std::invalid_argument("reduce_by_divide: numvec/denomvec does not represent an integer value.");
        }
    }
    
    // Check for equality between factors in numerator and denominator
    for (size_t i = 0; i < denomvec.size(); ++i) {
        int di = denomvec[i];
        auto it = std::find(numvec.begin(), numvec.end(), di);
        if (it != numvec.end()) {
            *it = 1;
            denomvec[i] = 1;
        }
    }
    
	std::erase_if(numvec, f);
	std::erase_if(denomvec, f);
    
    if (denomvec.empty()) {
        if (numvec.empty()) {
            numvec = {1};
            denomvec = {1};
            return;
        } else {
            denomvec = {1};
            return;
        }
    } else {
        if (numvec.empty()) {
            throw std::invalid_argument("numvec/denomvec does not represent an integer value.");
        }
    }
    
    // 20 first prime numbers
    std::vector<unsigned int> divisors = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71};
    
    for (size_t j = 0; j < denomvec.size(); ++j) {
        // Search decomposition of denomvec[j] as a product of prime numbers
        unsigned int denomj = denomvec[j];
        std::vector<unsigned int> denomj_decomposition;
        while (denomj > 1) {
            for (const unsigned int & d : divisors) {
                if (denomj % d == 0) {
                    denomj_decomposition.push_back(d);
                    denomj = denomj / d;
                    break;
                }
            }
        }
        std::sort(denomj_decomposition.begin(), denomj_decomposition.end());
        
        for (const unsigned int & d : denomj_decomposition) {
            // Find first element in numvec divisible by d
            for (size_t k = 0; k < numvec.size(); ++k) {
                if (numvec[k] % d == 0) {
                    numvec[k] = numvec[k] / d;
                    denomvec[j] = denomvec[j] / d;
                    break;
                }
            }
        }
    }
    
	std::erase_if(numvec, f);
	std::erase_if(denomvec, f);
    
    if (numvec.empty()) {
        if (!denomvec.empty()) {
            throw std::invalid_argument("numvec/denomvec does not represent an integer value.");
        }
        numvec = {1};
    }
    
    if (denomvec.empty()) {
        denomvec = {1};
    }
}


UINT64 multiplicity(const Irrep & alpha) {
	
	unsigned int n = alpha.n();
	unsigned int nrows = alpha.nrows();
    std::vector<unsigned int> alphaT = alpha.transpose();
	
	std::vector<unsigned int> arnum(n-1);
	for (size_t i=0; i<arnum.size(); ++i) {
        arnum[i] = i+2;
    }
    
    std::vector<unsigned int> denom;
    for (unsigned int i=0; i<nrows; ++i) {
        for (unsigned int j=0; j<alpha[i]; ++j) {
            int sum = alpha[i] + alphaT[j] - j - i - 1;
			denom.push_back(sum);
        }
    }
    
    reduce_by_divide(arnum, denom);
    
    if ((denom.size()>1) || (denom[0]>1)) {
		std::abort();
	}
	UINT64 falpha = 1;
	for (unsigned int i=0; i<arnum.size(); ++i) {
		falpha *= arnum[i];
	}
	return falpha;
}

} // namespace sun
