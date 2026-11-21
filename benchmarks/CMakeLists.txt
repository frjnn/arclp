# Table 1 of Yang (2018) on the Netlib problems.
add_executable(arclp_bench
    bench_main.c
    bench_report.c
    bench_solve.c
    paper_table.c)
target_link_libraries(arclp_bench PRIVATE arclp::arclp arclp::mat5
                                          arclp_compiler_options)

# Checks every solution against the original (un-presolved) problem.
add_executable(netlib_verify netlib_verify.c paper_table.c)
target_link_libraries(netlib_verify PRIVATE arclp::arclp arclp::mat5
                                            arclp_compiler_options)

# `cmake --build <dir> --target benchmark` runs the full table.
add_custom_target(benchmark
    COMMAND arclp_bench ${ARCLP_NETLIB_DIR} --repeat 3
            --csv ${CMAKE_BINARY_DIR}/netlib_results.csv
    DEPENDS arclp_bench
    USES_TERMINAL
    COMMENT "Running the Netlib benchmark")
