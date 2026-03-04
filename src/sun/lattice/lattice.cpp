// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include "lattice.h"

#include "../permutation/permutation.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <stdexcept>


namespace sun {

Lattice::Lattice(const std::string & filename)
{
	std::ifstream file;
	file.open(filename);

	std::vector<std::string> allLines;
	if (file.is_open()==true) {
		while (file.good())
		{
			std::string line;
			std::getline(file, line);
			if (!line.empty() && (line.back() == '\r')) {
				line.pop_back();
			}
			allLines.push_back(line);
		}
		file.close();
	} else {
		throw std::runtime_error("Problem opening lattice file. Aborting.");
	}

	if (allLines.size()==0) {
		throw std::runtime_error("Problem: lattice file is empty. Aborting.");
	}

	if (allLines[allLines.size()-1].size()==0) {
		allLines.pop_back();
	}
	
	/*
	std::cout << "------------------------------------------" << std::endl;
	std::cout << "Lattice file: " << filename << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	for (unsigned int i=0; i<allLines.size(); ++i) {
		std::cout << allLines[i] << std::endl;
	}
	std::cout << "------------------------------------------" << std::endl;
	*/
	
	unsigned int jsites = 0;
	unsigned int jbond = 0;
	unsigned int jcaching = 0;

	Ns_ = 0;
	nbonds_ = 0;
	ncaching_ = 0;
	
	for (unsigned int j=0; j<allLines.size(); ++j)
	{
		std::string line = allLines[j];
		std::string flag;

		flag = "[sites]=";
		if ((line.size()>=flag.size()+1) && (line.substr(0, flag.size())==flag)) {
			Ns_ = std::stoi(line.substr(flag.size()));
			jsites = j;
		}

		flag = "[interactions]=";
		if ((line.size()>=flag.size()+1) && (line.substr(0, flag.size())==flag)) {
			nbonds_ = std::stoi(line.substr(flag.size()));
			jbond = j;
		}
		
		flag = "[caching]=";
		if ((line.size()>=flag.size()+1) && (line.substr(0, flag.size())==flag)) {
			ncaching_ = std::stoi(line.substr(flag.size()));
			jcaching = j;
		}
	}

	if (jsites+Ns_+1!=jbond) {
		throw std::runtime_error("Problem: missing/too many site definition in lattice file.");
	}
	
	// new - with caching
	if (ncaching_==0) {
		// no caching request
		if ((jcaching>0) && (jbond+nbonds_!=jcaching)) {
			throw std::runtime_error("Problem: missing bonds in lattice file (A).");
		} else if ((jcaching==0) && (jbond+nbonds_>=allLines.size())) {
			throw std::runtime_error("Problem: missing bonds in lattice file (B).");
		}
	} else {
		if (jbond+nbonds_+1!=jcaching) {
			throw std::runtime_error("Problem: missing/too many bonds definition in lattice file.");
		}
		
		if (jcaching+ncaching_>=allLines.size()) {
			throw std::runtime_error("Problem: missing caching details in lattice file.");
		}
	}
	
	parse_sites(jsites, allLines);
	parse_bonds(jbond, allLines);
	//parse_caching_info(jcaching, allLines);
}


void Lattice::parse_sites(const unsigned int jsites, const std::vector<std::string>& allLines)
{
	sites_.resize(Ns_);
	
	for (unsigned int jj=0; jj<Ns_; ++jj) {
		
		const unsigned int j = jj + jsites + 1;
		std::string line = allLines[j];
		std::stringstream stream(line);
		sites_[jj].resize(3);
		int n;
		double x, y, z;
		stream >> n >> x >> y >> z;
		if (n!=static_cast<int>(jj)) {
			throw std::runtime_error(std::string("Problem: site n=") + std::to_string(n) + " not in order in lattice file.");
		}
		sites_[n][0] = x;
		sites_[n][1] = y;
		sites_[n][2] = z;
	}
}


void Lattice::parse_bonds(const unsigned int jbond, const std::vector<std::string>& allLines)
{
	bonds.resize(nbonds_);
	couplings_.resize(0);
	
	for (unsigned int jj=0; jj<nbonds_; ++jj)
	{
		const unsigned int j = jj + jbond + 1;
		std::string line = allLines[j];
		unsigned int cpt0 = 0; // counter for opening parenthesis
		unsigned int cpt1 = 0; // counter for closing parenthesis
		for (unsigned int k=0; k<line.size(); ++k) {
			if (line[k]=='(') {
				cpt0 += 1;
			} else if (line[k]==')') {
				cpt1 += 1;
			}
		}
		if (cpt0!=cpt1) {
			throw std::runtime_error(std::string("Problem: cycles not well defined at interaction ") + std::to_string(jj) + ". Parenthesis unmatch.");
		}

		// Start by reading coupling name
		if (line[0]=='(') {
			throw std::runtime_error(std::string("Problem: interaction ") + std::to_string(jj) + " should start with a coupling.");
		}
		unsigned int ck = 0;
		while (line[ck]!=' ') {
			ck += 1;
		}
		if (ck==0) {
			throw std::runtime_error(std::string("Problem: interaction ") + std::to_string(jj) + " should start with a coupling.");
		}
		std::string couplingName = line.substr(0, ck);
		
		bonds[jj].couplingName = couplingName;
		bonds[jj].couplingValue = 0.0;

		if (couplings_.size()==0) {
			couplings_.push_back(couplingName);
		} else {
			unsigned int c;
			for (c=0; c<couplings_.size(); ++c) {
				if (couplings_[c]==couplingName) {
					break;
				}
			}
			if (c==couplings_.size()) {
				couplings_.push_back(couplingName);
			}
		}

		std::vector<int> tmpcycle(0);
		unsigned int pos = ck+1;
		Permutation ptmp;

		if (line[pos]!='(') {
			throw std::runtime_error(std::string("Problem: interaction ") + std::to_string(jj) + ": syntax wrong.");
		}

		while (pos<line.size()) {
			if (line[pos]=='(') {
				tmpcycle.resize(0);
				pos += 1;
			} else if (line[pos]==',') {
				pos += 1;
			} else if (line[pos]==' ') {
				pos += 1;
			} else if (line[pos]==')') {
				// we have finished reading a cycle
				if (tmpcycle.size()==2) {
					ptmp.push_back(Transposition(tmpcycle));
				} else {
					throw std::runtime_error("Lattice constructor: here, we treat only transpositions, not p-cycles with p>2");
					ptmp.push_back(Cycle(tmpcycle));
				}
				pos += 1;
			} else {
				// read a site
				unsigned int le = 1;
				while ((line[pos+le]!=',') && (line[pos+le]!=')')) {
					le += 1;
				}
				unsigned int site = std::stoi(line.substr(pos, le));
				if (site>=Ns_) {
					throw std::runtime_error(std::string("Problem: interaction ") + std::to_string(jj) + ": site can not be >=Ns.");
				}
				tmpcycle.push_back(site);
				pos += le;
			}
		}
		
		// Reverse all cycles to satisfy the convention: right most cycle is
		// operated first, and should appear at the top of the vector
		std::reverse(ptmp.begin(), ptmp.end());
		bonds[jj].p = ptmp;
			
		bonds[jj].ops = toAdjacentTransposition(bonds[jj].p); // bonds[jj].ops is a std::vector<AdjacentTransposition>
		
		// id_ will order bonds by length, and then by starting index
		bonds[jj].id_ = (bonds[jj].ops[0].getk()) | (bonds[jj].ops.size()<<8);
	}
	
	std::sort(bonds.begin(), bonds.end(), [](const Bond& a, const Bond& b) {
		return a.id_ < b.id_;
	});
}


void Lattice::parse_caching_info(const unsigned int jcaching, const std::vector<std::string>& allLines)
{
	if (jcaching>0) {
		for (unsigned int jj=0; jj<ncaching_; ++jj) {
			const unsigned int j = jj + jcaching + 1;
			std::string line = allLines[j];
			std::stringstream stream(line);
			unsigned int bcid;
			stream >> bcid;
			if (bcid>=nbonds_) {
				throw std::runtime_error("caching bond index is not valid");
			}
			cache_indices_.push_back(bcid);
			bonds[bcid].cacheit_ = true;
		}
		
		for (size_t bond_index=0; bond_index<nbonds_; ++bond_index) {
			const unsigned int k = bonds[bond_index].ops[0].getk();
			
			for (size_t child_bond_index=bond_index+1; child_bond_index<nbonds_; ++child_bond_index) {
				const unsigned int kchild = bonds[child_bond_index].ops[0].getk();
				if (kchild==k) {
					bonds[bond_index].child_ = child_bond_index;
					bonds[child_bond_index].parent_ = bond_index;
					break;
				}
			}
		}
	}
}


void Lattice::print_sites() const
{
	std::cout << "------------------------------------------" << std::endl;
	std::cout << "Parsed sites: " << Ns_ << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	for (unsigned int j=0; j<Ns_; ++j) {
		std::cout << j << " (" << sites_[j][0] << ", " << sites_[j][1] << ", " << sites_[j][2] << ")" << std::endl;
	}
}

void Lattice::print_bonds() const
{
	std::cout << "------------------------------------------" << std::endl;
	std::cout << "Parsed bonds: " << nbonds_ << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	for (unsigned int j=0; j<nbonds_; ++j) {
		std::cout << std::right << std::setw(2) << j << ") [" 
				  << std::right << std::setw(6) << bonds[j].id_ << "] "
		          << bonds[j].couplingName << "=" 
				  << bonds[j].couplingValue << ": " 
				  << "[" << bonds[j].ops.size() << "] : " 
				  << bonds[j].p << " = " << bonds[j].ops << std::endl;
	}
}


void Lattice::print_bonds_light() const
{
	std::cout << "------------------------------------------" << std::endl;
	std::cout << "Parsed bonds: " << nbonds_ << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	for (unsigned int j=0; j<nbonds_; ++j) {
		std::cout << bonds[j].p << std::endl;
	}
}


void Lattice::print_bonds_light_with_decomposition() const
{
	std::cout << "------------------------------------------" << std::endl;
	std::cout << "Parsed bonds: " << nbonds_ << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	for (unsigned int j=0; j<nbonds_; ++j) {
		std::cout << std::right << std::setw(2) << j << ") {"
				  << (bonds[j].ops.size()+1)/2
				  << "} : " 
				  << bonds[j].p 
				  << " = "
				  << bonds[j].ops
				  << std::endl;
	}
}


void Lattice::print_caching_info() const
{
	std::cout << "------------------------------------------" << std::endl;
	std::cout << "Caching: " << ncaching_ << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	for (unsigned int j=0; j<ncaching_; ++j) {
		
		auto it = bonds[cache_indices_[j]].ops.begin();
		std::vector<AdjacentTransposition> sub(it, it + (bonds[cache_indices_[j]].ops.size()-1)/2+1);
		
		std::cout << std::right << std::setw(2) << j << ") " 
				  << "index=" << std::right << std::setw(2) << cache_indices_[j] << " : "
				  << "[" << bonds[j].ops.size() << "] : from " 
				  << bonds[cache_indices_[j]].p << " caching " << sub << std::endl;
	}
}


void Lattice::print_parents() const
{
	std::cout << "------------------------------------------" << std::endl;
	std::cout << "Parents: " << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	for (unsigned int j=0; j<nbonds_; ++j) {
		if (bonds[j].parent_!=-1) {
			std::cout << std::right << std::setw(2) << j << ") " 
					  << bonds[j].parent_ << std::endl;
		}
	}
}


void Lattice::print_children() const
{
	std::cout << "------------------------------------------" << std::endl;
	std::cout << "Caching lines: " << std::endl;
	std::cout << "------------------------------------------" << std::endl;
	for (unsigned int j=0; j<nbonds_; ++j) {
		const Bond* bond_ptr = &bonds[j];
		if ((bond_ptr->cacheit_==true) && (bond_ptr->parent_==-1) && (bond_ptr->child_>-1)) {
			// Follow the chain
			std::cout << j << std::flush;
			while (bond_ptr->child_>-1) {
				std::cout << " ---> " << bond_ptr->child_ << std::flush;
				bond_ptr = &bonds[bond_ptr->child_];
			}
			std::cout << std::endl;
		}
	}
}


void Lattice::assignCouplingValues(const std::map<std::string, double> & couplingsMap)
{
	for (unsigned int i=0; i<bonds.size(); ++i) {
		bonds[i].couplingValue = couplingsMap.at(bonds[i].couplingName);
	}
}


double Lattice::get_average_number_of_adjacent_transpositions_per_bond() const
{
	double avg = 0.0;
	for (unsigned int j=0; j<nbonds_; ++j) {
		avg += bonds[j].ops.size();
	}
	
	return avg/nbonds_;
}

} // namespace sun
