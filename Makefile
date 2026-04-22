
.PHONY: \
	all \
	build_mpi_eigvals build_mpi_eigvecs build_mpi_correlations \
	build_debug_mpi build_tests_mpi \
	clean

all: build_mpi_eigvals build_mpi_eigvecs build_mpi_correlations

#===================================================================
# 3 most important targets: MPI x {eigvals, eigvecs, correlations}
#===================================================================

build_mpi_eigvals:
	cmake -S . -B build/release_mpi_eigvals -DCMAKE_BUILD_TYPE=Release -DSUNED_USE_MPI=ON -DSUNED_EIGVALS=ON
	cmake --build build/release_mpi_eigvals
	ln -sf build/release_mpi_eigvals/main_mpi main_mpi_eigvals

build_mpi_eigvecs:
	cmake -S . -B build/release_mpi_eigvecs -DCMAKE_BUILD_TYPE=Release -DSUNED_USE_MPI=ON -DSUNED_EIGVECS=ON
	cmake --build build/release_mpi_eigvecs
	ln -sf build/release_mpi_eigvecs/main_mpi main_mpi_eigvecs

build_mpi_correlations:
	cmake -S . -B build/release_mpi_correlations -DCMAKE_BUILD_TYPE=Release -DSUNED_USE_MPI=ON -DSUNED_CORRELATIONS=ON
	cmake --build build/release_mpi_correlations
	ln -sf build/release_mpi_correlations/main_mpi main_mpi_correlations

#===================================================================
# Debug on MPI version
#===================================================================

build_debug_mpi:
	cmake -S . -B build/debug_mpi -DCMAKE_BUILD_TYPE=Debug -DSUNED_USE_MPI=ON
	cmake --build build/debug_mpi

#===================================================================
# Tests on MPI version
#===================================================================

build_tests_mpi:
	cmake -S . -B build/testing_mpi -DBUILD_TESTING=ON -DSUNED_USE_MPI=ON
	cmake --build build/testing_mpi

clean:
	rm -rf build
