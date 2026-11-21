# Installs the library, its public headers, a CMake package
# (find_package(arclp) -> arclp::arclp) and a pkg-config file.
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS arclp
    EXPORT arclpTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/arclp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

set(ARCLP_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/arclp)
install(EXPORT arclpTargets
    NAMESPACE arclp::
    DESTINATION ${ARCLP_CMAKE_DIR})
configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/arclpConfig.cmake.in
    ${PROJECT_BINARY_DIR}/arclpConfig.cmake
    INSTALL_DESTINATION ${ARCLP_CMAKE_DIR})
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/arclpConfigVersion.cmake
    COMPATIBILITY SameMajorVersion)
install(FILES
    ${PROJECT_BINARY_DIR}/arclpConfig.cmake
    ${PROJECT_BINARY_DIR}/arclpConfigVersion.cmake
    DESTINATION ${ARCLP_CMAKE_DIR})

if(ARCLP_BLAS_BACKEND STREQUAL "accelerate")
    set(ARCLP_PC_LIBS_PRIVATE "-framework Accelerate -lm")
elseif(ARCLP_BLAS_BACKEND STREQUAL "cblas")
    list(JOIN BLAS_LIBRARIES " " _arclp_blas_flags)
    set(ARCLP_PC_LIBS_PRIVATE "${_arclp_blas_flags} -lm")
else()
    set(ARCLP_PC_LIBS_PRIVATE "-lm")
endif()
configure_file(${PROJECT_SOURCE_DIR}/cmake/arclp.pc.in
    ${PROJECT_BINARY_DIR}/arclp.pc @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/arclp.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
