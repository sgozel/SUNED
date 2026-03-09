// Copyright 2026 Samuel GOZEL, GNU GPLv3

#ifndef SG_TIME_H
#define SG_TIME_H

#include <string>
#include <chrono>


void time(
	const std::chrono::time_point<std::chrono::high_resolution_clock>& tstart, 
	const std::chrono::time_point<std::chrono::high_resolution_clock>& tend,
	const std::string& def);


void time(
	const std::chrono::time_point<std::chrono::high_resolution_clock>& tstart, 
	const std::string& def);

#endif
