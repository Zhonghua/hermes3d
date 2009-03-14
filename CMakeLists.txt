project(h3d)

cmake_minimum_required(VERSION 2.6)

# default values

set(DEBUG       NO)
set(DEBUG_ORDER NO)
set(PROFILING   NO)

# Hermes3D related

# real/complex version of the library
set(REAL        YES)
set(COMPLEX     NO)
# supported element types
set(WITH_TETRA  YES)
set(WITH_HEX    YES)
set(WITH_PRISM  NO)

set(WITH_UMFPACK YES)
set(WITH_PETSC   NO)
set(WITH_PARDISO NO)

set(WITH_TESTS  YES)
set(WITH_TOOLS  NO)

# Doxygen related
set(DOXYGEN_BINARY doxygen)
set(DOXYGEN_CONFIG_FILE ${PROJECT_SOURCE_DIR}/doc/doxygen-config)

# allow to override default values via CMake.vars 
if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
    include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

# sanity checks

# need at least real or complex version
if(REAL)
    list(APPEND VERS "real")
endif(REAL)
if(COMPLEX)
    list(APPEND VERS "complex")
endif(COMPLEX)

list(LENGTH VERS VERS_LEN)
if(VERS_LEN GREATER 0)
    message(STATUS "Will build versions: ${VERS}")
else(VERS_LEN GREATER 0)
	message(FATAL_ERROR "No versions specified. Re-run with REAL and/or COMPLEX set to YES")
endif(VERS_LEN GREATER 0)

if(WITH_TETRA)
    list(APPEND ELEMTYPES "tetra")
endif(WITH_TETRA)
if(WITH_HEX)
    list(APPEND ELEMTYPES "hex")
endif(WITH_HEX)
if(WITH_PRISM)
    list(APPEND ELEMTYPES "prism")
endif(WITH_PRISM)

list(LENGTH ELEMTYPES ELEMTYPES_LEN)
if(ELEMTYPES_LEN GREATER 0)
    message(STATUS "Will build with element types: ${ELEMTYPES}")
else(ELEMTYPES_LEN GREATER 0)
    message(FATAL_ERROR "No element types selected. Re-run with WITH_TETRA, WITH_HEX or WITH_PRISM set to YES")
endif(ELEMTYPES_LEN GREATER 0)

# set internal library-wide variables
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/common)
set(MESHES_DIR ${PROJECT_SOURCE_DIR}/meshes)

# internals regarding build type
if(DEBUG)
	if(REAL)
	    set(HERMES_REAL_BIN hermes3d-debug)
	endif(REAL)
	if(COMPLEX) 
	    set(HERMES_CPLX_BIN hermes3d-cplx-debug)
	endif(COMPLEX)
	set(CMAKE_BUILD_TYPE Debug)
else(DEBUG)
	if(REAL)
	    set(HERMES_REAL_BIN hermes3d)
	endif(REAL)
	if(COMPLEX)
	    set(HERMES_CPLX_BIN hermes3d-cplx)
	endif(COMPLEX)
	set(CMAKE_BUILD_TYPE Release)
endif(DEBUG)

# If using PETSc or UMFPack, we need to enable Fortran support and look for BLAS and LAPACK
if(WITH_PETSC OR WITH_UMFPACK)
    enable_language(Fortran)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif(WITH_PETSC OR WITH_UMFPACK)

# find necessary packages
find_package(JUDY REQUIRED)


# linear solvers
if(WITH_PETSC)
	find_package(PETSC REQUIRED)
endif(WITH_PETSC)

if(WITH_UMFPACK)
    find_package(UMFPACK REQUIRED)
endif(WITH_UMFPACK)

if(WITH_PARDISO)
    find_package(PARDISO REQUIRED)
endif(WITH_PARDISO)


if(WITH_HDF5)
    find_package(HDF5 REQUIRED)
endif(WITH_HDF5)

if(WITH_MPI)
    find_package(MPI REQUIRED)
endif(WITH_MPI)

if(WITH_GLUT)
    find_package(GLUT REQUIRED)
endif(WITH_GLUT)

include_directories(${PROJECT_SOURCE_DIR})

#
add_subdirectory(src)
if(WITH_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif(WITH_TESTS)
if(WITH_TOOLS)
    add_subdirectory(tools)
endif(WITH_TOOLS)

# doc
add_subdirectory(doc)

add_custom_target(doc)
add_custom_command(
    SOURCE    ${DOXYGEN_CONFIG_FILE}
    COMMAND   ${DOXYGEN_BINARY}
    ARGS      ${DOXYGEN_CONFIG_FILE}
    TARGET    doc
    OUTPUTS   ${PROJECT_BINARY_DIR}/doc/html
)

add_custom_command(
    SOURCE    doc
    TARGET    doc
    DEPENDS   ${PROJECT_BINARY_DIR}/doc/html
)

add_custom_target(test-quick
    COMMAND   /usr/bin/ctest -E \"hnnd|qorder|lobatto|hex-h1-1face\"
)

# -- package --------

set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

set(CPACK_PACKAGE_VENDOR "hermes3d development team")
set(CPACK_SYSTEM_NAME "linux")

# source
set(CPACK_SOURCE_GENERATOR "TBZ2")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "hermes3d-${VERSION}")
set(CPACK_SOURCE_IGNORE_FILES 
    "${PROJECT_SOURCE_DIR}/.hg/"
    "${PROJECT_SOURCE_DIR}/.settings/"
    "${PROJECT_SOURCE_DIR}/.hgignore"
    "${PROJECT_SOURCE_DIR}/.project"
    "${PROJECT_SOURCE_DIR}/.cproject"
    "${PROJECT_SOURCE_DIR}/build/"
    "${PROJECT_SOURCE_DIR}/CMake.vars$"
    "${PROJECT_SOURCE_DIR}/utils/"
)

include(CPack)
