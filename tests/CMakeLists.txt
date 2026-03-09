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
# COMMON TEST SOURCES (ALWAYS BUILT)
#=======================================================================

set(COMMON_TEST_PROJECT_SOURCES
	${CMAKE_SOURCE_DIR}/src/common/time.cpp
	${CMAKE_SOURCE_DIR}/src/tmatrix/tmatrix.cpp
    ${CMAKE_SOURCE_DIR}/src/sun/utils/utils.cpp
)

set(COMMON_TEST_SOURCES
	test_datatypes.cpp
    test_falpha.cpp
    test_irrep.cpp
    test_tmatrix.cpp
    test_bsyt.cpp
)

add_executable(tests_common ${COMMON_TEST_SOURCES} ${COMMON_TEST_PROJECT_SOURCES})

add_dependencies(tests_common update_version)

target_include_directories(tests_common PRIVATE ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(tests_common PRIVATE
    GTest::gtest
    GTest::gtest_main
    BLAS::BLAS
    OpenMP::OpenMP_CXX
)

include(GoogleTest)
gtest_discover_tests(tests_common)

#=======================================================================
# MPI/NON-MPI DEDICATED TESTS
#=======================================================================

if(USE_MPI)

	#=======================================================================
	# MPI PROJECT SOURCES
	#=======================================================================
	
    set(MPI_PROJECT_SOURCES
		${COMMON_SOURCES}
        ${CMAKE_SOURCE_DIR}/src/ed_mpi/ed_solver_mpi.cpp
        ${CMAKE_SOURCE_DIR}/src/sun/heisenberg/hb_engine.cpp
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

    add_executable(tests_mpi ${MPI_TEST_SOURCES} ${MPI_PROJECT_SOURCES})

	add_dependencies(tests_mpi update_version)

    target_include_directories(tests_mpi PRIVATE 
								${CMAKE_SOURCE_DIR}/src
								${CMAKE_BINARY_DIR}/src/version)

    target_compile_definitions(tests_mpi PRIVATE SG_USE_MPI)

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

else()

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
	
	add_executable(tests_nonmpi ${NONMPI_TEST_SOURCES} ${NONMPI_PROJECT_SOURCES})
	
	add_dependencies(tests_nonmpi update_version)
	
	target_include_directories(tests_nonmpi PRIVATE 
									${CMAKE_SOURCE_DIR}/src
									${CMAKE_BINARY_DIR}/src/version)

	target_link_libraries(tests_nonmpi PRIVATE
		GTest::gtest
		GTest::gtest_main
		BLAS::BLAS
		OpenMP::OpenMP_CXX
	)

	include(GoogleTest)
	gtest_discover_tests(tests_nonmpi)
	
endif()
