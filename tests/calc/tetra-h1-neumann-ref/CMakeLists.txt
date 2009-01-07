project(calc-tetra-h1-neumann-ref)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	tetra-neumann.cc
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
	add_test(${PROJECT_NAME}-1 ${BIN} ${MESHES_DIR}/mesh3d/tetra1-ref.mesh3d 2)
else(FAST_TEST)
	foreach(a RANGE 2 10)
		add_test(${PROJECT_NAME}-1-${a} ${BIN} ${MESHES_DIR}/mesh3d/tetra1-ref.mesh3d ${a})
	endforeach(a)
endif(FAST_TEST)
