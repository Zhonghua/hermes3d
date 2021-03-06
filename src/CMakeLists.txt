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
	linsolver/pardiso.cc
	linsolver/petsc.cc
	linsolver/umfpack.cc
	${HERMES_COMMON_DIR}/callstack.cc
	${HERMES_COMMON_DIR}/trace.cc
	${HERMES_COMMON_DIR}/error.cc
	${HERMES_COMMON_DIR}/utils.cc
	${HERMES_COMMON_DIR}/timer.cc
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
	h3dconfig.h
	hermes3d.h
	integrals/h1.h
	matrix.h
	linsolver.h
	norm.h
	order.h
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
	${HERMES_COMMON_DIR}/arrayptr.h
	${HERMES_COMMON_DIR}/bitarray.h
	${HERMES_COMMON_DIR}/map.h
	${HERMES_COMMON_DIR}/maphs.h
	${HERMES_COMMON_DIR}/mapord.h
	${HERMES_COMMON_DIR}/utils.h
	${HERMES_COMMON_DIR}/error.h
	${HERMES_COMMON_DIR}/trace.h
	${HERMES_COMMON_DIR}/timer.h
	${HERMES_COMMON_DIR}/callstack.h
)

set(INC_SPACE
	space/h1.h
	space/hcurl.h
)

set(INC_SHAPESET
	shapeset/common.h
	shapeset/h1lobattotetra.h
	shapeset/h1lobattohex.h
	shapeset/hcurllobattohex.h
	shapeset/hex.h
	shapeset/tetra.h
	shapeset/refmapss.h
	shapeset/lobatto.h
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
	output/graph.h
)

set(INC_INTEGRALS
	integrals/h1.h
	integrals/hcurl.h
)


macro(LIBRARY_SETTINGS BIN)
	if(WITH_PETSC)
		include_directories(${PETSC_INCLUDE_DIRS})
		target_link_libraries(${BIN} ${PETSC_LIBRARIES})
		target_link_libraries(${BIN} ${PETSC_LIBS})

       	target_link_libraries(${BIN} ${BLAS_LIBRARIES})
   	    target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
   	endif(WITH_PETSC)

	if(WITH_UMFPACK)
		include_directories(${UMFPACK_INCLUDE_DIR} ${AMD_INCLUDE_DIR})
		target_link_libraries(${BIN} ${UMFPACK_LIBRARY} ${AMD_LIBRARY})

       	target_link_libraries(${BIN} ${BLAS_LIBRARIES})
   	    target_link_libraries(${BIN} ${LAPACK_LIBRARIES})
	endif(WITH_UMFPACK)

	if(WITH_HDF5)
	endif(WITH_HDF5)
endmacro(LIBRARY_SETTINGS)


# real version
if(REAL) 
    add_library(${HERMES_REAL_BIN} SHARED ${SRC})
	LIBRARY_SETTINGS(${HERMES_REAL_BIN})
endif(REAL)

# complex version
if(COMPLEX)
    add_library(${HERMES_CPLX_BIN} SHARED ${SRC})
    set_target_properties(${HERMES_CPLX_BIN} PROPERTIES COMPILE_FLAGS -DCOMPLEX)
	LIBRARY_SETTINGS(${HERMES_CPLX_BIN})
endif(COMPLEX)

#
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/h3dconfig.h
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
install(FILES ${INC_INTEGRALS} DESTINATION include/${PROJECT_NAME}/integrals)

# install library
if(REAL)
    install(TARGETS ${HERMES_REAL_BIN} LIBRARY DESTINATION lib/${PROJECT_NAME})
endif(REAL)

if(COMPLEX)
    install(TARGETS ${HERMES_CPLX_BIN} LIBRARY DESTINATION lib/${PROJECT_NAME})
endif(COMPLEX)

