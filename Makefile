
.PHONY: all \
        build build_numa \
        build_mpi build_mpi_numa \
        build_debug build_debug_mpi 
        build_test build_tests_numa \
        build_test_mpi build_tests_mpi_numa \
        build_tests_all \
        clean

all: build

build:
	cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
	cmake --build build/release
	ln -sf build/release/main main

build_numa:
	cmake -S . -B build/release_numa -DCMAKE_BUILD_TYPE=Release -DUSE_NUMA=ON
	cmake --build build/release_numa
	ln -sf build/release_numa/main main_numa

build_mpi:
	cmake -S . -B build/release_mpi -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON
	cmake --build build/release_mpi
	ln -sf build/release_mpi/main_mpi main_mpi

build_mpi_numa:
	cmake -S . -B build/release_mpi_numa -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON -DUSE_NUMA=ON
	cmake --build build/release_mpi_numa
	ln -sf build/release_mpi_numa/main_mpi main_mpi_numa

build_debug:
	cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug
	cmake --build build/debug

build_debug_mpi:
	cmake -S . -B build/debug_mpi -DCMAKE_BUILD_TYPE=Debug -DUSE_MPI=ON
	cmake --build build/debug_mpi

build_tests:
	cmake -S . -B build/testing -DBUILD_TESTING=ON
	cmake --build build/testing

build_tests_numa:
	cmake -S . -B build/testing_numa -DBUILD_TESTING=ON -DUSE_NUMA=ON
	cmake --build build/testing_numa

build_tests_mpi:
	cmake -S . -B build/testing_mpi -DBUILD_TESTING=ON -DUSE_MPI=ON
	cmake --build build/testing_mpi

build_tests_mpi_numa:
	cmake -S . -B build/testing_mpi_numa -DBUILD_TESTING=ON -DUSE_MPI=ON -DUSE_NUMA=ON
	cmake --build build/testing_mpi_numa

build_tests_all: build_tests build_tests_numa build_tests_mpi build_tests_mpi_numa

clean:
	rm -rf build

