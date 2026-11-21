# Optional MATLAB v5 .mat reader used by the tests and the benchmark.
find_package(ZLIB REQUIRED)

add_library(arclp_mat5 STATIC
    src/mat5_array.c
    src/mat5_element.c
    src/mat5_file.c
    src/mat5_lp.c)
add_library(arclp::mat5 ALIAS arclp_mat5)
target_include_directories(arclp_mat5
    PUBLIC include
    PRIVATE src)
target_link_libraries(arclp_mat5
    PUBLIC ZLIB::ZLIB
    PRIVATE arclp_compiler_options)
