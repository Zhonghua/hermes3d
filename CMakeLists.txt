project(h3d)

cmake_minimum_required(VERSION 2.4)

include(CMake.vars)
include_directories(${PROJECT_SOURCE_DIR})

# need either real or complex version
if(NOT REAL AND NOT COMPLEX)
	message(FATAL_ERROR "Specify either REAL and/or COMPLEX version of Hermes3D in global CMake.vars")
endif(NOT REAL AND NOT COMPLEX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(HERMES_COMMON_DIR ${PROJECT_SOURCE_DIR}/common)
set(MESHES_DIR ${PROJECT_SOURCE_DIR}/meshes)

# build type
if(DEBUG)
	set(PETSC_ARCH ${PETSC_ARCH_DEBUG})
	if(REAL)
		set(HERMES_BIN hermes3d-debug)
	endif(REAL)
	if(COMPLEX) 
		set(HERMES_CPLX_BIN hermes3d-cplx-debug)
	endif(COMPLEX)
	set(CMAKE_BUILD_TYPE Debug)
else(DEBUG)
	set(PETSC_ARCH ${PETSC_ARCH_OPT})
	if(REAL)
		set(HERMES_BIN hermes3d)
	endif(REAL)
	if(COMPLEX) 
		set(HERMES_CPLX_BIN hermes3d-cplx)
	endif(COMPLEX)
	set(CMAKE_BUILD_TYPE Release)
endif(DEBUG)

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

