build:
	cmake -S . -B build/release -DCMAKE_BUILD_TYPE=Release
	cmake --build build/release

build_mpi:
	cmake -S . -B build/mpi -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=ON
	cmake --build build/mpi

debug:
	cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug
	cmake --build build/debug

debug_mpi:
	cmake -S . -B build/mpi_debug -DCMAKE_BUILD_TYPE=Debug -DUSE_MPI=ON
	cmake --build build/mpi_debug

clean:
	rm -rf build

test:
	cmake -S . -B build/testing -DBUILD_TESTING=ON
	cmake --build build/testing

test_mpi:
	cmake -S . -B build/testing_mpi -DBUILD_TESTING=ON -DUSE_MPI=ON
	cmake --build build/testing_mpi
