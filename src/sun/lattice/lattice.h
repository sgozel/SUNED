// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SUN_LATTICE_H
#define SUN_LATTICE_H

#include <string>
#include <vector>
#include <map>

#include "bond.h"

namespace sun {

class Lattice
{
public:
	Lattice(const std::string & filename);
	Lattice() = default;

	unsigned int get_Ns() const {return Ns_;};
	unsigned int get_nbonds() const {return nbonds_;};
	std::vector<std::string> get_couplingNames() const {return couplings_;};
	void assignCouplingValues(const std::map<std::string, double> &);
	
	void print_sites() const;
	void print_bonds() const;
	void print_bonds_light() const;
	void print_bonds_light_with_decomposition() const;
	void print_caching_info() const;
	void print_parents() const;
	void print_children() const;
	
	double get_average_number_of_adjacent_transpositions_per_bond() const;

private:
	void parse_sites(const unsigned int jsites, const std::vector<std::string>& allLines);
	void parse_bonds(const unsigned int jbond, const std::vector<std::string>& allLines);
	void parse_caching_info(const unsigned int jcaching, const std::vector<std::string>& allLines);
	
public:
	std::vector<Bond> bonds;

private:
	unsigned int Ns_;
	unsigned int nbonds_;
	std::vector<std::vector<double> > sites_; // (x, y, z) coordinates of each site
	std::vector<std::string> couplings_;
	
	// caching
	unsigned int ncaching_;
	std::vector<unsigned int> cache_indices_;
};

} // namespace sun

#endif
