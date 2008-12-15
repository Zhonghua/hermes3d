project(tetra-h1-dirichlet)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	tetra-dirichlet.cc
	${HERMES_COMMON_DIR}/timer.cc
	${HERMES_COMMON_DIR}/error.cc
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 

include_directories(${hermes3d_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

#
# Tests
#

enable_testing()

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
if(FAST_TEST)
	add_test(tetra-h1-dirichlet-1 ${BIN} ${MESHES_DIR}/mesh3d/tetra1.mesh3d 4)
	add_test(tetra-h1-dirichlet-2 ${BIN} ${MESHES_DIR}/mesh3d/tetra2.mesh3d 3)
	add_test(tetra-h1-dirichlet-3 ${BIN} ${MESHES_DIR}/mesh3d/tetra4.mesh3d 2)
	add_test(tetra-h1-dirichlet-4 ${BIN} ${MESHES_DIR}/mesh3d/tetra8.mesh3d 2)
else(FAST_TEST)
	foreach(a RANGE 4 10)
		add_test(tetra-h1-dirichlet-1-${a} ${BIN} ${MESHES_DIR}/mesh3d/tetra1.mesh3d ${a})
	endforeach(a)

	foreach(a RANGE 3 10)
		add_test(tetra-h1-dirichlet-2-${a} ${BIN} ${MESHES_DIR}/mesh3d/tetra2.mesh3d ${a})
	endforeach(a)
	
	foreach(a RANGE 2 10)
		add_test(tetra-h1-dirichlet-3-${a} ${BIN} ${MESHES_DIR}/mesh3d/tetra4.mesh3d ${a})
	endforeach(a)
	
	foreach(a RANGE 2 10)
		add_test(tetra-h1-dirichlet-4-${a} ${BIN} ${MESHES_DIR}/mesh3d/tetra8.mesh3d ${a})
	endforeach(a)
endif(FAST_TEST)
