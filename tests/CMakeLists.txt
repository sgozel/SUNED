#=======================================================================
# IMPORT LATTICE FILES
#=======================================================================

configure_file(../latticefiles/HB_chain_12_OBC.lattice HB_chain_12_OBC.lattice)
configure_file(../latticefiles/HB_chain_12_PBC.lattice HB_chain_12_PBC.lattice)
configure_file(../latticefiles/HB_chain_15_OBC.lattice HB_chain_15_OBC.lattice)
configure_file(../latticefiles/HB_chain_15_PBC.lattice HB_chain_15_PBC.lattice)

#=======================================================================
# IMPORT TEST DATA FILES
#=======================================================================

configure_file(testdata/TEST_DATA_MULTIPLICITY.json TEST_DATA_MULTIPLICITY.json)
configure_file(testdata/TEST_DATA_ENERGY_SU3.json TEST_DATA_ENERGY_SU3.json)

#=======================================================================
# MPI/NON-MPI DEDICATED TESTS
#=======================================================================

include(GoogleTest)

if(SUNED_USE_MPI)

	#=======================================================================
	# MPI PROJECT SOURCES
	#=======================================================================
	
    set(MPI_PROJECT_SOURCES
		${COMMON_SOURCES}
        ${CMAKE_SOURCE_DIR}/src/ed_mpi/ed_solver_mpi.cpp
        ${CMAKE_SOURCE_DIR}/src/sun/heisenberg_mpi/hb_fund_engine_mpi.cpp
        ${CMAKE_SOURCE_DIR}/src/sun/heisenberg_mpi/hb_fund_matrix_engine_mpi.cpp
    )

	#=======================================================================
	# MPI TEST SOURCES
	#=======================================================================

    set(MPI_TEST_SOURCES
		main_mpi_test.cpp
        test_energy_mpi.cpp
    )
	
	#=======================================================================
	# MPI TEST EXECUTABLE
	#=======================================================================
	
    add_executable(tests_mpi ${MPI_TEST_SOURCES} ${MPI_PROJECT_SOURCES})

	add_dependencies(tests_mpi update_version)

    target_include_directories(tests_mpi PRIVATE 
							   ${CMAKE_SOURCE_DIR}/src
							   ${CMAKE_BINARY_DIR}/src/version
							   ${CMAKE_CURRENT_SOURCE_DIR})

    target_compile_definitions(tests_mpi PRIVATE SG_USE_MPI)
    if(SUNED_USE_NUMA)
		target_compile_definitions(tests_mpi PRIVATE SG_USE_NUMA)
	endif()

    target_link_libraries(tests_mpi PRIVATE
        GTest::gtest
        BLAS::BLAS
        OpenMP::OpenMP_CXX
        MPI::MPI_CXX
    )
    
    set(MPI_TEST_NP 2 CACHE STRING "Number of MPI processes for tests")
	add_test(
		NAME tests_mpi
		COMMAND mpirun --oversubscribe -np ${MPI_TEST_NP} $<TARGET_FILE:tests_mpi>
	)
	set_tests_properties(tests_mpi PROPERTIES
		ENVIRONMENT "OMP_NUM_THREADS=1"
	)
	
	# do not use gtest_discover_tests, which shadows the add_test with the mpirun command
	#gtest_discover_tests(tests_mpi)
	
else()
	
	#=======================================================================
	# COMMON TEST PROJECT SOURCES
	#=======================================================================
	
	set(COMMON_TEST_PROJECT_SOURCES
		${CMAKE_SOURCE_DIR}/src/common/time.cpp
		${CMAKE_SOURCE_DIR}/src/tmatrix/tmatrix.cpp
		${CMAKE_SOURCE_DIR}/src/sun/utils/utils.cpp
		${CMAKE_SOURCE_DIR}/src/sun/syt_usage/bsyt_usage.cpp
	)
	
	#=======================================================================
	# COMMON TEST SOURCES
	#=======================================================================
	
	set(COMMON_TEST_SOURCES
		test_datatypes.cpp
		test_falpha.cpp
		test_irrep.cpp
		test_tmatrix.cpp
		test_bsyt.cpp
		test_bsyt_usage.cpp
		test_vsyt.cpp
		test_vsyt_usage.cpp
	)
	
	#=======================================================================
	# COMMON TEST EXECUTABLE
	#=======================================================================
	
	add_executable(tests_common ${COMMON_TEST_SOURCES} ${COMMON_TEST_PROJECT_SOURCES})

	add_dependencies(tests_common update_version)

	target_include_directories(tests_common PRIVATE 
							   ${CMAKE_SOURCE_DIR}/src
							   ${CMAKE_BINARY_DIR}/src/version
							   ${CMAKE_CURRENT_SOURCE_DIR})
	
	if(SUNED_USE_NUMA)
		target_compile_definitions(tests_common PRIVATE SG_USE_NUMA)
	endif()

	target_link_libraries(tests_common PRIVATE
		GTest::gtest
		GTest::gtest_main
		BLAS::BLAS
		OpenMP::OpenMP_CXX
	)

	gtest_discover_tests(tests_common)
	
	#=======================================================================
	# NON-MPI PROJECT SOURCES
	#=======================================================================

	set(NONMPI_PROJECT_SOURCES
		${COMMON_SOURCES}
		${CMAKE_SOURCE_DIR}/src/ed/ed_solver.cpp
		${CMAKE_SOURCE_DIR}/src/sun/heisenberg/hb_fund_engine.cpp
		${CMAKE_SOURCE_DIR}/src/sun/heisenberg/hb_fund_matrix_engine.cpp
		${CMAKE_SOURCE_DIR}/src/sun/heisenberg/hb_fund_mf_engine.cpp
	)

	#=======================================================================
	# NON-MPI TEST SOURCES
	#=======================================================================

	set(NONMPI_TEST_SOURCES
		test_energy.cpp
	)
	
	#=======================================================================
	# NON-MPI TEST EXECUTABLE
	#=======================================================================
	
	add_executable(tests_nonmpi ${NONMPI_TEST_SOURCES} ${NONMPI_PROJECT_SOURCES})
	
	add_dependencies(tests_nonmpi update_version)
	
	target_include_directories(tests_nonmpi PRIVATE 
							   ${CMAKE_SOURCE_DIR}/src
							   ${CMAKE_BINARY_DIR}/src/version
							   ${CMAKE_CURRENT_SOURCE_DIR})
	
	if(SUNED_USE_NUMA)
		target_compile_definitions(tests_nonmpi PRIVATE SG_USE_NUMA)
	endif()
	
	target_link_libraries(tests_nonmpi PRIVATE
		GTest::gtest
		GTest::gtest_main
		BLAS::BLAS
		OpenMP::OpenMP_CXX
	)
	
	gtest_discover_tests(tests_nonmpi)
	
endif()
