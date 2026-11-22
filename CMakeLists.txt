cmake_minimum_required(VERSION 3.21)

project(arclp
    VERSION 1.0.0
    DESCRIPTION "Arc-search infeasible interior-point solver for linear programs"
    LANGUAGES C)

# ---- Options ---------------------------------------------------------------
option(ARCLP_BUILD_TESTS "Build the test suite" ${PROJECT_IS_TOP_LEVEL})
option(ARCLP_BUILD_BENCHMARKS "Build the Netlib benchmark" ${PROJECT_IS_TOP_LEVEL})
option(ARCLP_BUILD_EXAMPLES "Build the examples" ${PROJECT_IS_TOP_LEVEL})
option(ARCLP_NATIVE "Optimize for the build machine (-march=native)" ON)
option(ARCLP_WERROR "Treat compiler warnings as errors" OFF)
set(ARCLP_SANITIZE "" CACHE STRING
    "Comma-separated sanitizers, e.g. address,undefined (empty: none)")
set(ARCLP_BLAS "auto" CACHE STRING
    "Dense kernels: auto, accelerate, cblas or none (portable C)")
set_property(CACHE ARCLP_BLAS PROPERTY STRINGS auto accelerate cblas none)
set(ARCLP_METIS "auto" CACHE STRING
    "Also try METIS nested dissection: auto (if found), on or off")
set_property(CACHE ARCLP_METIS PROPERTY STRINGS auto on off)
set(ARCLP_OPENMP "auto" CACHE STRING
    "Parallel loops with OpenMP: auto (if found), on or off")
set_property(CACHE ARCLP_OPENMP PROPERTY STRINGS auto on off)
set(ARCLP_NETLIB_DIR "${PROJECT_SOURCE_DIR}/data/arcLP" CACHE PATH
    "Directory holding the lp_*.mat Netlib files (scripts/fetch_netlib.sh)")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(GNUInstallDirs)
include(ArclpCompilerOptions)
include(ArclpBlas)
include(ArclpMetis)
include(ArclpOpenMP)

# ---- Targets ---------------------------------------------------------------
add_subdirectory(src)

if(ARCLP_BUILD_TESTS OR ARCLP_BUILD_BENCHMARKS)
    add_subdirectory(extras/mat5)
endif()

if(ARCLP_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(ARCLP_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(ARCLP_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

include(ArclpInstall)
