# Build-only interface targets carrying warnings, optimization and
# sanitizer flags.  They are linked PRIVATE and never exported.
include(CheckCCompilerFlag)

add_library(arclp_compiler_options INTERFACE)

if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(arclp_compiler_options INTERFACE
        -Wall -Wextra -Wpedantic -Wshadow -Wstrict-prototypes
        -Wmissing-prototypes -Wcast-qual -Wvla
        $<$<BOOL:${ARCLP_WERROR}>:-Werror>)
elseif(MSVC)
    target_compile_options(arclp_compiler_options INTERFACE
        /W4 $<$<BOOL:${ARCLP_WERROR}>:/WX>)
endif()

if(ARCLP_NATIVE)
    check_c_compiler_flag(-march=native ARCLP_HAS_MARCH_NATIVE)
    if(ARCLP_HAS_MARCH_NATIVE)
        target_compile_options(arclp_compiler_options INTERFACE -march=native)
    endif()
endif()

if(ARCLP_SANITIZE)
    target_compile_options(arclp_compiler_options INTERFACE
        -fsanitize=${ARCLP_SANITIZE} -fno-omit-frame-pointer
        -fno-sanitize-recover=all)
    target_link_options(arclp_compiler_options INTERFACE
        -fsanitize=${ARCLP_SANITIZE})
endif()

# glibc hides POSIX clocks under a strict -std=c11.
if(NOT APPLE AND NOT WIN32)
    target_compile_definitions(arclp_compiler_options INTERFACE
        _POSIX_C_SOURCE=200809L)
endif()
