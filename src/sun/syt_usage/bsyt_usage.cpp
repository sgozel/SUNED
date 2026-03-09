// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "bsyt_usage.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <utility>

#include "../utils/utils.h"

namespace sun {

std::vector<tbSYT> get_SYT(const Irrep & alpha) {
	
	unsigned int n = alpha.n();
	unsigned int nl = alpha.nrows();
	std::vector<unsigned int> alphaT = alpha.transpose();
	UINT64 falpha = multiplicity(alpha);
	
	std::vector<tbSYT> Y(falpha);
	
	// construct first SYT
    tbSYT y0;
    unsigned int cptel = 0;
    for (unsigned int j=0; j<alpha[0]; ++j) {
        for (unsigned int row=0; row<alphaT[j]; ++row) {
            y0.set(cptel, row);
            cptel += 1;
        }
    }
	
	Y[0] = y0;
    UINT64 s = 1;
    
    while (s<falpha) {
        tbSYT y(Y[s-1]);

        std::vector<unsigned int> lbd(nl+1, 0);
        lbd[0] = 1;

        unsigned int j = 1;
        while ((j<n) && (y[j]>=y[j-1])) {
            lbd[y[j]] += 1;
            j += 1;
        }
        lbd[y[j]] += 1;

        unsigned int t = lbd[y[j]+1];

        unsigned int i = nl;
        while (lbd[i-1] != t) {
            i -= 1;
        }
        y.set(j, i-1);
        lbd[i-1] -= 1;

        t = j;
        unsigned int k = 0;
        while (k<t) {
            unsigned int r = 1;
            while (lbd[r-1]>0) {
                y.set(k, r-1);
                lbd[r-1] -= 1;
                k += 1;
                r += 1;
            }
        }
        Y[s] = y;
        s += 1;
    }
    
    return Y;
}

std::vector<int> get_column(const tbSYT syt, const Irrep& alpha) {
	
	const unsigned int n = alpha.n();
	const unsigned int nr = alpha.nrows();
	
	std::vector<int> temp(nr, 0);
	temp[0] = 1;
	std::vector<int> cy(n, 0);
	
	for (unsigned int i=1; i<n; ++i) {
		int rowi = syt.get(i);
		cy[i] = temp[rowi];
		temp[rowi] += 1;
	}
	return cy;
}


std::pair<int, int> get_column_k_k_plus_one(const tbSYT syt, const unsigned int k) {
	std::pair<int, int> cols;
	cols.first = 0;
	cols.second = 0;
	const int rowk = syt.get(k);
	const int rowkk = syt.get(k+1);
	for (unsigned int i=0; i<k; ++i) {
		const int rowi = syt.get(i);
		if (rowi==rowk) {
			cols.first += 1;
		}
		if (rowi==rowkk) {
			cols.second += 1;
		}
	}
	return cols;
}


int get_axial_distance(const tbSYT y, const std::vector<int> & cy, const int i, const int j) {
	return cy[i] - y.get(i) - cy[j] + y.get(j);
}


void fullsimplify_development(std::vector<tbSYT> & ydev,
                              std::vector<double> & coeffdev,
                              const UINT64 sorted_up_to)
{
    if (ydev.size() <= 1) {
		return;
	}
    
    std::vector<tbSYT> c(ydev);
    
    // sort the unsorted portion (from sorted_up_to onwards)
    if (sorted_up_to < c.size()) {
        std::sort(c.begin() + sorted_up_to, c.end());
        // if there is a sorted prefix, merge it with the newly sorted suffix
        if (sorted_up_to > 0) {
            std::inplace_merge(c.begin(), c.begin() + sorted_up_to, c.end());
        }
    }
    
    // extract unique values
    auto last = std::unique(c.begin(), c.end());
    c.erase(last, c.end());
    
    // sum coefficients of equal SYTs
    std::vector<double> coeff_sum(c.size(), 0.0);
    for (size_t i = 0; i < ydev.size(); ++i) {
        auto it = std::lower_bound(c.begin(), c.end(), ydev[i]);
        coeff_sum[it - c.begin()] += coeffdev[i];
    }
    
    // remove elements with abs(coeff) < prec
    std::vector<tbSYT>  ydev_out;
    std::vector<double> coeff_out;
    ydev_out.reserve(c.size());
    coeff_out.reserve(c.size());
    const double prec = 1.0e-13;
    for (size_t i = 0; i < c.size(); ++i) {
        if (std::abs(coeff_sum[i]) >= prec) {
            ydev_out.push_back(c[i]);
            coeff_out.push_back(coeff_sum[i]);
        }
    }
    ydev = std::move(ydev_out);
    coeffdev = std::move(coeff_out);
}


void develop_consecutive_number(const Irrep & alpha,
								const std::vector<tbSYT> & ydev_in,
                                const std::vector<double> & coeffdev_in,
                                std::vector<tbSYT> & ydev_out,
                                std::vector<double> & coeffdev_out,
                                const int & k)
{
    UINT64 Ny = ydev_in.size();

    ydev_out = ydev_in;
    coeffdev_out = coeffdev_in;

    ydev_out.resize(2*Ny);
    coeffdev_out.resize(2*Ny);

    UINT64 count = Ny;

    for (UINT64 i=0; i<Ny; ++i) {
        
        const tbSYT& y = ydev_in[i];
        
        const int rowk = y.get(k);
        const int rowkk = y.get(k+1);
        
        if (rowk!=rowkk) {	
			std::pair<int, int> cy = get_column_k_k_plus_one(y, k);
			if (cy.first==cy.second) {
                coeffdev_out[i] *= (-1.0);
            } else {
                double rho = 1.0/static_cast<double>( rowk + cy.second - cy.first - rowkk ); // we have introduced the minus sign therein
                ydev_out[count] = y;
                ydev_out[count].exchange(k, k+1);
                coeffdev_out[count] = coeffdev_out[i] * std::sqrt(1.0-rho*rho);
                coeffdev_out[i] *= rho; // no minus sign, because rho is inverse axial distance from k+1 to k
                count += 1;
            }
        }
    }
    ydev_out.resize(count);
    coeffdev_out.resize(count);
}


void develop_consecutive_number_inplace(const Irrep & alpha,
										std::vector<tbSYT> & ydev,
										std::vector<double> & coeffdev,
										const int & k)
{
    UINT64 Ny = ydev.size();
    ydev.reserve(2*Ny);
    coeffdev.reserve(2*Ny);
    UINT64 count = Ny;

    for (UINT64 i=0; i<Ny; ++i) {
        
        const tbSYT& y = ydev[i];
        const int rowk = y.get(k);
        const int rowkk = y.get(k+1);
        
        if (rowk!=rowkk) {	
			const std::pair<int, int> cy = get_column_k_k_plus_one(y, k);
			if (cy.first==cy.second) {
                coeffdev[i] *= (-1.0);
            } else {
                double rho = 1.0/static_cast<double>( rowk + cy.second - cy.first - rowkk ); // we have introduced the minus sign therein
                ydev.push_back(y);
                ydev[count].exchange(k, k+1);
                coeffdev.push_back(coeffdev[i] * std::sqrt(1.0-rho*rho));
                coeffdev[i] *= rho; // no minus sign, because rho is inverse axial distance from k+1 to k
                count += 1;
            }
        }
    }
}

} // namespace sun
