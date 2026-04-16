cmake_minimum_required(VERSION 3.21)
project(SUNED VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g3")
set(COMMON_COMPILE_OPTIONS -Wall)
# for debug mode, invoke cmake as:
# cmake ../src/ -DCMAKE_BUILD_TYPE=Debug

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)

#=======================================================================
# OPTIONS
#=======================================================================

option(SUNED_EIGVALS "Build eigenvalue extraction" ON)
option(SUNED_EIGVECS "Build eigenvector extraction" OFF)
option(SUNED_CORRELATIONS "Build correlations extraction" OFF)
option(SUNED_USE_MPI "Build the MPI target (main_mpi) instead of the default target (main)" ON)
option(SUNED_USE_NUMA "Enable NUMA-aware allocation" OFF)
option(SUNED_LANCZOS_TWO_VECTORS "Build ED solver using two Lanczos vectors" OFF)
option(SUNED_USE_VSYT "Use basic vector-SYT implementation" OFF)
option(SUNED_STORE_COLUMNS "Store column SYTs" OFF)
set(SUNED_BSYT_NBITS "2" CACHE STRING "Number of bits per box (1, 2, 3, 4)")
set_property(CACHE SUNED_BSYT_NBITS PROPERTY STRINGS "1" "2" "3" "4")
set(SUNED_BSYT_CONTAINER_SIZE "8" CACHE STRING "Number of bytes per SYT (2, 4, 8)")
set_property(CACHE SUNED_BSYT_CONTAINER_SIZE PROPERTY STRINGS "2" "4" "8")

message(STATUS "SUNED_USE_VSYT            = ${SUNED_USE_VSYT}")

if(NOT SUNED_USE_VSYT)
    if(NOT SUNED_BSYT_NBITS MATCHES "^[1-4]$")
        message(FATAL_ERROR "SUNED_BSYT_NBITS should be 1, 2, 3 or 4 for SU(N) with 2<=N<=16. Got: ${SUNED_BSYT_NBITS}.")
    endif()
    if(NOT SUNED_BSYT_CONTAINER_SIZE MATCHES "^(2|4|8)$")
        message(FATAL_ERROR "SUNED_BSYT_CONTAINER_SIZE should be 2, 4 or 8. Got: ${SUNED_BSYT_CONTAINER_SIZE}.")
    endif()
    message(STATUS "SUNED_BSYT_NBITS          = ${SUNED_BSYT_NBITS}")
    message(STATUS "SUNED_BSYT_CONTAINER_SIZE = ${SUNED_BSYT_CONTAINER_SIZE}")
endif()

if(SUNED_EIGVECS)
    set(SUNED_EIGVALS OFF)
endif()

if(SUNED_CORRELATIONS)
    set(SUNED_EIGVALS OFF)
    set(SUNED_EIGVECS ON)
endif()

message(STATUS "SUNED_EIGVALS             = ${SUNED_EIGVALS}")
message(STATUS "SUNED_EIGVECS             = ${SUNED_EIGVECS}")
message(STATUS "SUNED_CORRELATIONS        = ${SUNED_CORRELATIONS}")
message(STATUS "SUNED_USE_MPI             = ${SUNED_USE_MPI}")
message(STATUS "SUNED_USE_NUMA            = ${SUNED_USE_NUMA}")

#=======================================================================
# SOURCE GIT INFO
#=======================================================================

execute_process(
    COMMAND git rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND git log -1 --format=%h
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
    COMMAND git log -1 --format=%ci
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_TIMESTAMP
    OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "CMAKE_SOURCE_DIR : ${CMAKE_SOURCE_DIR}")
message(STATUS "Git branch: ${GIT_BRANCH}")
message(STATUS "Git hash: ${GIT_COMMIT_HASH}")
message(STATUS "Git timestamp: ${GIT_COMMIT_TIMESTAMP}")

add_custom_target(update_version ALL
    COMMAND ${CMAKE_COMMAND}
        -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
        -DBINARY_DIR=${CMAKE_BINARY_DIR}
        -P ${CMAKE_SOURCE_DIR}/cmake/UpdateVersion.cmake
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

#=======================================================================
# CONFIGURE FILES
#=======================================================================

configure_file(data.json data.json)
include_directories(${CMAKE_BINARY_DIR}/src/version)

#=======================================================================
# MAIN TARGET
#=======================================================================

set(COMMON_SOURCES
    ${CMAKE_SOURCE_DIR}/src/common/time.cpp
    ${CMAKE_SOURCE_DIR}/src/tmatrix/tmatrix.cpp
    ${CMAKE_SOURCE_DIR}/src/sun/permutation/permutation.cpp
    ${CMAKE_SOURCE_DIR}/src/sun/lattice/lattice.cpp
    ${CMAKE_SOURCE_DIR}/src/sun/utils/utils.cpp
    ${CMAKE_SOURCE_DIR}/src/sun/heisenberg/hb_engine.cpp
)

if(NOT SUNED_USE_VSYT)
    list(APPEND COMMON_SOURCES ${CMAKE_SOURCE_DIR}/src/sun/syt_usage/bsyt_usage.cpp)
endif()

if(SUNED_USE_MPI)
    find_package(MPI REQUIRED)

    set(MPI_SOURCES
        ${CMAKE_SOURCE_DIR}/src/main_mpi.cpp
        ${CMAKE_SOURCE_DIR}/src/ed_mpi/ed_solver_mpi.cpp
        ${CMAKE_SOURCE_DIR}/src/sun/heisenberg_mpi/hb_fund_engine_mpi.cpp
        ${CMAKE_SOURCE_DIR}/src/sun/heisenberg_mpi/hb_fund_matrix_engine_mpi.cpp
    )

    add_executable(main_mpi ${MPI_SOURCES} ${COMMON_SOURCES})
    add_dependencies(main_mpi update_version)

    target_compile_options(main_mpi PRIVATE ${COMMON_COMPILE_OPTIONS})
    target_compile_definitions(main_mpi PRIVATE SG_USE_MPI)
    if(SUNED_EIGVALS)
        target_compile_definitions(main_mpi PRIVATE SG_EIGVALS)
    endif()
    if(SUNED_USE_NUMA)
		target_compile_definitions(main_mpi PRIVATE SG_USE_NUMA)
	endif()
    if(SUNED_USE_VSYT)
        target_compile_definitions(main_mpi PRIVATE SG_USE_VSYT)
    else()
        target_compile_definitions(main_mpi PRIVATE SG_BSYT_NBITS=${SUNED_BSYT_NBITS})
        target_compile_definitions(main_mpi PRIVATE SG_BSYT_CONTAINER_SIZE=${SUNED_BSYT_CONTAINER_SIZE})
    endif()
    if(SUNED_STORE_COLUMNS)
        target_compile_definitions(main_mpi PRIVATE SG_STORE_COLUMNS)
    endif()

    target_include_directories(main_mpi PRIVATE 
							   ${CMAKE_SOURCE_DIR}/src/
							   ${CMAKE_BINARY_DIR}/src/version)
    
    target_link_libraries(main_mpi PRIVATE 
				BLAS::BLAS 
				OpenMP::OpenMP_CXX 
				MPI::MPI_CXX)

else()
    set(MAIN_SOURCES
        ${CMAKE_SOURCE_DIR}/src/main.cpp
        ${CMAKE_SOURCE_DIR}/src/ed/ed_solver.cpp
        ${CMAKE_SOURCE_DIR}/src/sun/heisenberg/hb_fund_engine.cpp
        ${CMAKE_SOURCE_DIR}/src/sun/heisenberg/hb_fund_mf_engine.cpp
        ${CMAKE_SOURCE_DIR}/src/sun/heisenberg/hb_fund_matrix_engine.cpp
    )

    add_executable(main ${MAIN_SOURCES} ${COMMON_SOURCES})
    add_dependencies(main update_version)

    target_compile_options(main PRIVATE ${COMMON_COMPILE_OPTIONS})
    if(SUNED_EIGVALS)
        target_compile_definitions(main PRIVATE SG_EIGVALS)
    endif()
    if(SUNED_EIGVECS)
        target_compile_definitions(main PRIVATE SG_EIGVECS)
    endif()
    if(SUNED_LANCZOS_TWO_VECTORS)
        target_compile_definitions(main PRIVATE SG_LANCZOS_TWO_VECTORS)
    endif()
    if(SUNED_USE_NUMA)
		target_compile_definitions(main PRIVATE SG_USE_NUMA)
	endif()
    if(SUNED_USE_VSYT)
        target_compile_definitions(main PRIVATE SG_USE_VSYT)
    else()
        target_compile_definitions(main PRIVATE SG_BSYT_NBITS=${SUNED_BSYT_NBITS})
        target_compile_definitions(main PRIVATE SG_BSYT_CONTAINER_SIZE=${SUNED_BSYT_CONTAINER_SIZE})
    endif()
    if(SUNED_STORE_COLUMNS)
        target_compile_definitions(main PRIVATE SG_STORE_COLUMNS)
    endif()

    target_include_directories(main PRIVATE 
							   ${CMAKE_SOURCE_DIR}/src/
							   ${CMAKE_BINARY_DIR}/src/version)
    
    target_link_libraries(main PRIVATE 
			BLAS::BLAS 
			OpenMP::OpenMP_CXX)

endif()

#=======================================================================
# MULTIPLICITY TARGET
#=======================================================================

add_executable(main_falpha
    src/main_falpha.cpp
    src/sun/utils/utils.cpp
)
add_dependencies(main_falpha update_version)

target_compile_options(main_falpha PRIVATE ${COMMON_COMPILE_OPTIONS})
target_include_directories(main_falpha PRIVATE 
						   ${CMAKE_SOURCE_DIR}/src/
						   ${CMAKE_BINARY_DIR}/src/version)

#=======================================================================
# LATTICE TARGET
#=======================================================================

add_executable(main_lattice
    src/main_lattice.cpp
    src/sun/lattice/lattice.cpp
    src/sun/permutation/permutation.cpp
)
add_dependencies(main_lattice update_version)

target_compile_options(main_lattice PRIVATE -Wall)
target_include_directories(main_lattice PRIVATE 
						   ${CMAKE_SOURCE_DIR}/src/
						   ${CMAKE_BINARY_DIR}/src/version)

#=======================================================================
# SYTS TARGET
#=======================================================================

set(SYT_TARGET_SOURCES
    ${CMAKE_SOURCE_DIR}/src/main_syts.cpp
    ${CMAKE_SOURCE_DIR}/src/sun/permutation/permutation.cpp
    ${CMAKE_SOURCE_DIR}/src/sun/utils/utils.cpp
    ${CMAKE_SOURCE_DIR}/src/common/time.cpp
)

if(NOT SUNED_USE_VSYT)
    list(APPEND SYT_TARGET_SOURCES ${CMAKE_SOURCE_DIR}/src/sun/syt_usage/bsyt_usage.cpp)
endif()

add_executable(main_syts
    ${SYT_TARGET_SOURCES}
)
add_dependencies(main_syts update_version)

target_compile_options(main_syts PRIVATE ${COMMON_COMPILE_OPTIONS})
if(SUNED_USE_VSYT)
    target_compile_definitions(main_syts PRIVATE SG_USE_VSYT)
else()
    target_compile_definitions(main_syts PRIVATE SG_BSYT_NBITS=${SUNED_BSYT_NBITS})
    target_compile_definitions(main_syts PRIVATE SG_BSYT_CONTAINER_SIZE=${SUNED_BSYT_CONTAINER_SIZE})
endif()
if(SUNED_USE_NUMA)
	target_compile_definitions(main_syts PRIVATE SG_USE_NUMA)
endif()

target_include_directories(main_syts PRIVATE 
						   ${CMAKE_SOURCE_DIR}/src/
						   ${CMAKE_BINARY_DIR}/src/version)

target_link_libraries(main_syts PRIVATE OpenMP::OpenMP_CXX)

#=======================================================================
# TESTS
#=======================================================================

if(PROJECT_IS_TOP_LEVEL)
	option(BUILD_TESTING "Build tests" OFF)  # override default
    include(CTest)
    if(BUILD_TESTING)
        message(STATUS "Building tests")
        find_package(GTest REQUIRED)
        add_subdirectory(tests)
    endif()
endif()
