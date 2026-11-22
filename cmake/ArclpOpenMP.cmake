# Optional OpenMP loop parallelism.  Sets ARCLP_USE_OPENMP and
# ARCLP_OPENMP_LIBRARIES (link paths, safe to export).

set(ARCLP_USE_OPENMP OFF)
set(ARCLP_OPENMP_LIBRARIES "")
if(NOT ARCLP_OPENMP STREQUAL "off")
    if(CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
        # Apple's clang has no OpenMP runtime; use Homebrew's libomp.
        find_path(ARCLP_LIBOMP_INCLUDE_DIR omp.h
            HINTS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include)
        find_library(ARCLP_LIBOMP_LIBRARY omp
            HINTS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
        if(ARCLP_LIBOMP_INCLUDE_DIR AND ARCLP_LIBOMP_LIBRARY)
            set(OpenMP_C_FLAGS
                "-Xpreprocessor -fopenmp -isystem ${ARCLP_LIBOMP_INCLUDE_DIR}")
            set(OpenMP_C_LIB_NAMES omp)
            set(OpenMP_omp_LIBRARY ${ARCLP_LIBOMP_LIBRARY})
        endif()
    endif()
    find_package(OpenMP COMPONENTS C)
    if(OpenMP_C_FOUND)
        set(ARCLP_USE_OPENMP ON)
        set(ARCLP_OPENMP_LIBRARIES ${OpenMP_C_LIBRARIES})
    elseif(ARCLP_OPENMP STREQUAL "on")
        message(FATAL_ERROR "ARCLP_OPENMP=on but OpenMP was not found")
    endif()
endif()
message(STATUS "arclp OpenMP: ${ARCLP_USE_OPENMP}")
