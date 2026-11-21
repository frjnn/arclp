# The solver.  Sources are compiled once into an object library, from which
# the public library (static or shared, per BUILD_SHARED_LIBS) and a static
# library exposing the internal headers for white-box tests are linked.

add_library(arclp_objects OBJECT
    api/arclp.c
    api/solve_reduced.c
    chol/amd_degree.c
    chol/amd_element.c
    chol/amd_init.c
    chol/amd_order.c
    chol/amd_supervar.c
    chol/assemble.c
    chol/chol_solver.c
    chol/etree.c
    chol/graph.c
    chol/numeric.c
    chol/solve.c
    chol/supernodes.c
    chol/symbolic.c
    chol/symbolic_assembly.c
    chol/symbolic_order.c
    chol/symbolic_structure.c
    dense/dense_cholesky.c
    ${ARCLP_DENSE_SOURCE}
    ipm/arc_angle.c
    ipm/arc_bounds.c
    ipm/arc_component.c
    ipm/arc_directions.c
    ipm/arc_rule.c
    ipm/ipm.c
    ipm/ipm_run.c
    ipm/mehrotra_directions.c
    ipm/mehrotra_length.c
    ipm/mehrotra_rule.c
    ipm/start_point.c
    presolve/presolve.c
    presolve/ps_matrix.c
    presolve/reduced.c
    presolve/rule_implied_free.c
    presolve/rules_basic.c
    presolve/substitutions.c
    sparse/csc.c
    util/arena.c)

set_target_properties(arclp_objects PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    C_VISIBILITY_PRESET hidden)
target_include_directories(arclp_objects
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${ARCLP_BLAS_INCLUDE_DIRS})
target_compile_definitions(arclp_objects PRIVATE ARCLP_BUILDING
    $<$<BOOL:${BUILD_SHARED_LIBS}>:ARCLP_SHARED>)
target_link_libraries(arclp_objects PRIVATE arclp_compiler_options)

find_library(ARCLP_LIBM m)
set(_arclp_link_deps ${ARCLP_BLAS_LIBRARIES})
if(ARCLP_LIBM)
    list(APPEND _arclp_link_deps ${ARCLP_LIBM})
endif()

# ---- Public library --------------------------------------------------------
add_library(arclp $<TARGET_OBJECTS:arclp_objects>)
add_library(arclp::arclp ALIAS arclp)
target_include_directories(arclp PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_link_libraries(arclp PRIVATE ${_arclp_link_deps}
                            $<BUILD_INTERFACE:arclp_compiler_options>)
if(BUILD_SHARED_LIBS)
    target_compile_definitions(arclp INTERFACE ARCLP_SHARED)
endif()
set_target_properties(arclp PROPERTIES
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR})

# ---- Internal library for white-box tests ----------------------------------
if(ARCLP_BUILD_TESTS)
    add_library(arclp_internal STATIC $<TARGET_OBJECTS:arclp_objects>)
    target_include_directories(arclp_internal PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(arclp_internal PUBLIC ${_arclp_link_deps})
endif()
