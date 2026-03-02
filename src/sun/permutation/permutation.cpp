// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "permutation.h"

#include <algorithm>
#include <stdexcept>


namespace sun {

Cycle::Cycle()
{
    cycle_.resize(0);
}


Cycle::Cycle(const std::initializer_list<int> & init)
: cycle_{init}
{
    if (cycle_.size()==1) {
        throw std::runtime_error("1-cycles are identities. Remove them.");
    }
    for (unsigned int i=1; i<cycle_.size(); ++i) {
        for (unsigned int j=0; j<i; ++j) {
            if (cycle_[i]==cycle_[j]) {
                throw std::runtime_error("Problem: element appears twice in cycle.");
            }
        }
    }
}


Cycle::Cycle(const std::vector<int> & init)
{
    cycle_ = init;
    if (cycle_.size()==1) {
        throw std::runtime_error("1-cycles are identities. Remove them.");
    }
    for (unsigned int i=1; i<cycle_.size(); ++i) {
        for (unsigned int j=0; j<i; ++j) {
            if (cycle_[i]==cycle_[j]) {
                throw std::runtime_error("Problem: element appears twice in cycle.");
            }
        }
    }
}


std::vector<Transposition> Cycle::toTransposition() const
{
    unsigned int n = cycle_.size();
    std::vector<Transposition> transpos(0);
    for (unsigned int i=0; i<n-1; ++i) {
        transpos.push_back(Transposition(cycle_[0], cycle_[i+1]));
    }
    return transpos;
}


std::vector<AdjacentTransposition> Cycle::toAdjacentTransposition() const
{
    std::vector<AdjacentTransposition> adjaTranspos(0);
	
    std::vector<Transposition> transpovec = this->toTransposition();
	
    for (unsigned int i=0; i<transpovec.size(); ++i) {
        std::vector<AdjacentTransposition> adjatranspoveci = transpovec[i].toAdjacentTransposition();
        for (unsigned int j=0; j<adjatranspoveci.size(); ++j) {
            adjaTranspos.push_back(adjatranspoveci[j]);
        }
    }
    return adjaTranspos;
}


Transposition::Transposition()
{
    cycle_.resize(0);
}


Transposition::Transposition(const int i, const int j)
: Cycle({std::min(i, j), std::max(i, j)})
{
}

Transposition::Transposition(const std::initializer_list<int> & init)
: Cycle(init)
{
	if (cycle_.size()!=2) {
        throw std::runtime_error("Wrong declaration of Transposition.");
    }
	if (cycle_[0]==cycle_[1]) {
		throw std::runtime_error("Problem: element appears twice in Transposition.");
	}
    if (cycle_[0]>cycle_[1]) {
		std::swap(cycle_[0], cycle_[1]);
	}
}

Transposition::Transposition(const std::vector<int> & init)
: Cycle(init)
{
	if (cycle_.size()!=2) {
        throw std::runtime_error("Wrong declaration of Transposition.");
    }
	if (cycle_[0]==cycle_[1]) {
		throw std::runtime_error("Problem: element appears twice in Transposition.");
	}
    if (cycle_[0]>cycle_[1]) {
		std::swap(cycle_[0], cycle_[1]);
	}
}


std::vector<AdjacentTransposition> Transposition::toAdjacentTransposition() const
{
    unsigned int minP = std::min(cycle_[0], cycle_[1]);
    unsigned int maxP = std::max(cycle_[0], cycle_[1]);
    unsigned int nbt = 2*(maxP - minP) - 1;

    std::vector<AdjacentTransposition> adjaTranspos(0);
    if (nbt==1) {
		adjaTranspos.push_back(AdjacentTransposition(minP));
	} else {
		for (unsigned int t=minP; t<maxP; ++t) {
			adjaTranspos.push_back(AdjacentTransposition(t));
		}
		for (unsigned int t=0; t<(nbt-1)/2; ++t) {
			adjaTranspos.push_back(AdjacentTransposition(maxP-2-t));
		}
	}
    return adjaTranspos;
}


AdjacentTransposition::AdjacentTransposition()
{
    cycle_.resize(0);
}


AdjacentTransposition::AdjacentTransposition(const int k)
: Transposition(k, k+1)
{
}


std::vector<AdjacentTransposition> AdjacentTransposition::toAdjacentTransposition() const
{
    return std::vector<AdjacentTransposition>{*this};
}


std::vector<AdjacentTransposition> toAdjacentTransposition(const Permutation& p)
{
    std::vector<AdjacentTransposition> adjaTranspos;

    for (unsigned int i=0; i<p.size(); ++i) {
        std::vector<AdjacentTransposition> ati = p[i].toAdjacentTransposition();
        for (unsigned int j=0; j<ati.size(); ++j) {
            adjaTranspos.push_back(ati[j]);
        }
    }
    return adjaTranspos;
}


std::ostream& operator<<(std::ostream& os, const Cycle& cycle)
{
    os << "(";
    for (unsigned int l=0; l<cycle.cycle_.size()-1; ++l) {
        os << cycle.cycle_[l] << ", ";
    }
    os << cycle.cycle_[cycle.cycle_.size()-1] << ")";
    return os;
}


void print(const Permutation & perm)
{
	for (unsigned int t=0; t<perm.size(); ++t) {
		std::cout << perm[t] << std::flush;
	}
	std::cout << std::endl;
}


void print(const std::vector<AdjacentTransposition> & adjaTranspos)
{
	for (unsigned int t=0; t<adjaTranspos.size(); ++t) {
		std::cout << adjaTranspos[t] << std::flush;
	}
	std::cout << std::endl;
}

} // namespace sun
