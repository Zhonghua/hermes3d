project(hermes3d)

include(CMake.vars)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})		# for config.h

# source files for the hermes3d library
set(SRC
	discretization.cc
	filter.cc
	function.cc
	mesh.cc
	loader/mesh3d.cc
	loader/hdf5.cc
	matrix.cc
	norm.cc
	output/gmsh.cc
	output/vtk.cc
	output/graph.cc
	precalc.cc
	qsort.cc
	quadstd.cc
	refdomain.cc
	refmap.cc
	shapeset.cc
	shapeset/lobatto.cc
	shapeset/h1lobattotetra.cc
	shapeset/h1lobattotetradx.cc
	shapeset/h1lobattotetrady.cc
	shapeset/h1lobattotetradz.cc
	shapeset/h1lobattohex.cc
	shapeset/hcurllobattohex.cc
	shapeset/refmapss.cc
	solution.cc
	space.cc
	space/h1.cc
	space/hcurl.cc
	transform.cc
	traverse.cc
	linsolver.cc
	linsolver/pardiso.cc
	linsolver/petsc.cc
	linsolver/umfpack.cc
	${HERMES_COMMON_DIR}/trace.cc
	${HERMES_COMMON_DIR}/error.cc
	${HERMES_COMMON_DIR}/utils.cc
)

# include files that will be installed
set(INC
	asmlist.h
	common.h
	discretization.h
	filter.h
	function.h
	mesh.h
	meshloader.h
	discretization.h
	hermes3d.h
	integrals/h1.h
	matrix.h
	linsolver.h
	norm.h
	output.h
	precalc.h
	quad.h
	quadstd.h
	refdomain.h
	refmap.h
	shapeset.h
	solution.h
	space.h
	transform.h
	traverse.h
)

set(INC_COMMON
	${HERMES_COMMON_DIR}/array.h
	${HERMES_COMMON_DIR}/bitarray.h
	${HERMES_COMMON_DIR}/map.h
	${HERMES_COMMON_DIR}/maphs.h
	${HERMES_COMMON_DIR}/mapord.h
	${HERMES_COMMON_DIR}/utils.h
)

set(INC_SPACE
	space/h1.h
	space/hcurl.h
)

set(INC_SHAPESET
	shapeset/common.h
	shapeset/h1lobattotetra.h
	shapeset/h1lobattohex.h
)

set(INC_LOADER
	loader/mesh3d.h
	loader/hdf5.h
)

set(INC_LINSOLVER
	linsolver/pardiso.h
	linsolver/petsc.h
	linsolver/umfpack.h
)

set(INC_OUTPUT
	output/gmsh.h
	output/vtk.h
)


# real version
if(REAL) 
    add_library(${HERMES_REAL_BIN} STATIC ${SRC})
endif(REAL)

# complex version
if(COMPLEX)
    add_library(${HERMES_CPLX_BIN} STATIC ${SRC})
    set_target_properties(${HERMES_CPLX_BIN} PROPERTIES COMPILE_FLAGS -DCOMPLEX)
endif(COMPLEX)

#
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)


#
# INSTALL
#

# install header files
install(FILES ${INC}           DESTINATION include/${PROJECT_NAME})
install(FILES ${INC_COMMON}    DESTINATION include/${PROJECT_NAME}/common)
install(FILES ${INC_SPACE}     DESTINATION include/${PROJECT_NAME}/space)
install(FILES ${INC_SHAPESET}  DESTINATION include/${PROJECT_NAME}/shapeset)
install(FILES ${INC_LOADER}    DESTINATION include/${PROJECT_NAME}/loader)
install(FILES ${INC_LINSOLVER} DESTINATION include/${PROJECT_NAME}/linsolver)
install(FILES ${INC_OUTPUT}    DESTINATION include/${PROJECT_NAME}/output)
install(FILES ${INC_ADAPT}     DESTINATION include/${PROJECT_NAME}/adapt)

# install library
if(REAL)
    install(TARGETS ${HERMES_REAL_BIN} ARCHIVE DESTINATION lib/${PROJECT_NAME})
endif(REAL)

if(COMPLEX)
    install(TARGETS ${HERMES_CPLX_BIN} ARCHIVE DESTINATION lib/${PROJECT_NAME})
endif(COMPLEX)

