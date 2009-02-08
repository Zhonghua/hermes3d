project(h3d)

cmake_minimum_required(VERSION 2.4)

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

set(USE_BLAS    YES)
set(USE_LAPACK  YES)
set(USE_FORTRAN gfortran)

set(USE_UMFPACK YES)
set(USE_PETSC   NO)
set(USE_PARDISO NO)

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
if(NOT REAL AND NOT COMPLEX)
	message(FATAL_ERROR "Specify either REAL and/or COMPLEX version of Hermes3D in global CMake.vars")
endif(NOT REAL AND NOT COMPLEX)

if(NOT WITH_TETRA AND NOT WITH_HEX AND NOT WITH_PRISM)
    message(FATAL_ERROR "No support of any type of elements")
endif(NOT WITH_TETRA AND NOT WITH_HEX AND NOT WITH_PRISM)

# set internal library-wide variables
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/common)
set(MESHES_DIR ${PROJECT_SOURCE_DIR}/meshes)

# internals regarding build type
if(DEBUG)
	set(PETSC_ARCH ${PETSC_ARCH_DEBUG})
	if(REAL)
	    set(HERMES_REAL_BIN hermes3d-debug)
	endif(REAL)
	if(COMPLEX) 
	    set(HERMES_CPLX_BIN hermes3d-cplx-debug)
	endif(COMPLEX)
	set(CMAKE_BUILD_TYPE Debug)
else(DEBUG)
	set(PETSC_ARCH ${PETSC_ARCH_OPT})
	if(REAL)
	    set(HERMES_REAL_BIN hermes3d)
	endif(REAL)
	if(COMPLEX)
	    set(HERMES_CPLX_BIN hermes3d-cplx)
	endif(COMPLEX)
	set(CMAKE_BUILD_TYPE Release)
endif(DEBUG)


# find necessary packages
find_package(JUDY REQUIRED)

# linear solvers
if(USE_PETSC)
   find_package(PETSC REQUIRED)
endif(USE_PETSC)

if(USE_UMFPACK)
    find_package(UMFPACK REQUIRED)
endif(USE_UMFPACK)

if(USE_PARDISO)
    find_package(PARDISO REQUIRED)
endif(USE_PARDISO)

if(USE_HDF5)
    find_package(HDF5 REQUIRED)
	find_package(ZLIB REQUIRED)
endif(USE_HDF5)

if(USE_MPI)
    find_package(MPI)
endif(USE_MPI)

if(USE_GLUT)
    find_package(GLUT)
endif(USE_GLUT)

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
    COMMAND   /usr/bin/ctest -E \"hnnd|qorder|lobatto\"
)
