// Copyright 2026 Samuel GOZEL, GNU GPLv3

#include <iostream>

template<class SYT>
void print_syt(const SYT& syt, const Irrep & alpha, const bool & zeroBased) {
	unsigned int n = alpha.n();
	unsigned int nr = alpha.nrows();
	std::vector<std::vector<unsigned int>> boxes(nr, std::vector<unsigned int>(0));
	for (unsigned int i=0; i<n; ++i) {
		int ri = syt.get(i);
		boxes[ri].push_back(i);
	}
	
	int shift = (zeroBased ? 0 : 1);
	
	for (unsigned int i=0; i<nr; ++i) {
		for (size_t j=0; j<boxes[i].size()-1; ++j) {
			std::cout << boxes[i][j]+shift << " & " << std::flush;
		}
		std::cout << boxes[i][boxes[i].size()-1]+shift << " \\\\" << std::endl;
	}
}
