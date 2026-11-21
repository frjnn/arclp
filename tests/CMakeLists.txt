# Unit tests reach into the internals (arclp_internal); integration tests
# use only the public API (arclp::arclp).  Netlib cases read their data
# from ARCLP_NETLIB_DIR and are skipped when it is absent.

function(arclp_add_test name source library label)
    add_executable(${name} ${source})
    target_include_directories(${name} PRIVATE helpers)
    target_link_libraries(${name} PRIVATE ${library} arclp::mat5
                                          arclp_compiler_options)
    add_test(NAME ${name} COMMAND ${name})
    set_tests_properties(${name} PROPERTIES
        LABELS ${label}
        ENVIRONMENT "NETLIB_DIR=${ARCLP_NETLIB_DIR}")
endfunction()

arclp_add_test(test_arc_angle unit/test_arc_angle.c arclp_internal unit)
arclp_add_test(test_dense unit/test_dense.c arclp_internal unit)
arclp_add_test(test_ipm_driver unit/test_ipm_driver.c arclp_internal unit)
arclp_add_test(test_cholesky unit/test_cholesky.c arclp_internal unit)
arclp_add_test(test_presolve unit/test_presolve.c arclp_internal unit)
arclp_add_test(test_solver integration/test_solver.c arclp::arclp integration)

find_package(Threads REQUIRED)
arclp_add_test(test_concurrency integration/test_concurrency.c arclp::arclp
               integration)
target_link_libraries(test_concurrency PRIVATE Threads::Threads)
if(ARCLP_BLAS_BACKEND STREQUAL "accelerate")
    target_compile_definitions(test_concurrency PRIVATE ARCLP_TEST_ACCELERATE)
    target_link_libraries(test_concurrency PRIVATE ${ARCLP_ACCELERATE})
endif()
