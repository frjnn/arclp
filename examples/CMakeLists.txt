add_executable(solve_small solve_small.c)
target_link_libraries(solve_small PRIVATE arclp::arclp arclp_compiler_options)
