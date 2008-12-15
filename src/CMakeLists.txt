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

macro(ADD_H3D_LIBRARY BIN)
	add_library(${BIN} STATIC ${SRC})

	#
	# linear solvers
	#

	# PETSc
	if(USE_PETSC)
		find_package(PETSC REQUIRED)
		if(PETSC_FOUND)
			include_directories(${PETSC_INCLUDE_DIRS})
			target_link_libraries(${BIN} ${PETSC_LIBRARIES})
		endif(PETSC_FOUND)
	endif(USE_PETSC)
	
	# UMFPACK
	if(USE_UMFPACK)
		find_package(UMFPACK REQUIRED)
		if(UMFPACK_FOUND)
			include_directories(
				${UMFPACK_INCLUDE_DIR}
				${AMD_INCLUDE_DIR}
			)
			target_link_libraries(${BIN}
				${UMFPACK_LIBRARY}
				${AMD_LIBRARY}
			)
		endif(UMFPACK_FOUND)
	endif(USE_UMFPACK)
	
	# PARDISO
	if(USE_PARDISO)
		find_package(PARDISO REQUIRED)
		if(PARDISO_FOUND)
			target_link_libraries(${BIN} ${PARDISO_LIBRARY})
		endif(PARDISO_FOUND)
	endif(USE_PARDISO)
	
	# Judy
	find_package(JUDY REQUIRED)
	if(JUDY_FOUND)
		include_directories(${JUDY_INCLUDE_DIR})
		target_link_libraries(${BIN} ${JUDY_LIBRARY})
	endif(JUDY_FOUND)
	
	# HDF5
	if(USE_HDF5)
		find_package(HDF5 REQUIRED)
		if(HDF5_FOUND)
			include_directories(${HDF5_INCLUDE_DIR})
			target_link_libraries(${BIN} ${HDF5_LIBRARY})
		endif(HDF5_FOUND)
	endif(USE_HDF5)
	
	# MPI
	if(USE_MPI)
		find_package(MPI)
		if(MPI_FOUND)
			include_directories(${MPI_INCLUDE_DIR})
			target_link_libraries(${BIN} ${MPI_LIBRARY})
		endif(MPI_FOUND)
	endif(USE_MPI)
	
	# GLUT
	if(USE_GLUT)
		find_package(GLUT)
		if(GLUT_FOUND)
			include_directories(${GLUT_INCLUDE_DIR})
			target_link_libraries(${BIN} ${GLUT_LIBRARY})
		endif(GLUT_FOUND)
	endif(USE_GLUT)
endmacro(ADD_H3D_LIBRARY)

# real version
if(REAL)
	ADD_H3D_LIBRARY(${HERMES_BIN})
endif(REAL)
# complex version
if(COMPLEX)
	ADD_H3D_LIBRARY(${HERMES_CPLX_BIN})
	set_target_properties(${HERMES_CPLX_BIN}
		PROPERTIES
		COMPILE_FLAGS -DCOMPLEX)
endif(COMPLEX)

# need at least one linear solver
if(NOT USE_PETSC AND NOT USE_UMFPACK AND NOT USE_PARDISO)
	message(FATAL_ERROR "You need to specify at least one linear solver library in global CMake.vars")
endif(NOT USE_PETSC AND NOT USE_UMFPACK AND NOT USE_PARDISO)

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
install(
	TARGETS ${HERMES_BIN}
	ARCHIVE DESTINATION lib/${PROJECT_NAME}
)
endif(REAL)

if(COMPLEX)
install(
	TARGETS ${HERMES_CPLX_BIN}
	ARCHIVE DESTINATION lib/${PROJECT_NAME}
)
endif(COMPLEX)

