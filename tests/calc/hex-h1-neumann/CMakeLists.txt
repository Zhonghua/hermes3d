project(hex-h1-neumann)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	hex-neumann.cc
	${HERMES_COMMON_DIR}/timer.cc
	${HERMES_COMMON_DIR}/error.cc
	${HERMES_COMMON_DIR}/utils.cc
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

add_test(hex-h1-neumann-1    ${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d 2 2 2)
add_test(hex-h1-neumann-sc-1 ${BIN} ${MESHES_DIR}/mesh3d/hex1-sc.mesh3d 2 2 2)
add_test(hex-h1-neumann-2x   ${BIN} ${MESHES_DIR}/mesh3d/hex2x.mesh3d 2 2 2)
add_test(hex-h1-neumann-2y   ${BIN} ${MESHES_DIR}/mesh3d/hex2y.mesh3d 2 2 2)
add_test(hex-h1-neumann-2z   ${BIN} ${MESHES_DIR}/mesh3d/hex2z.mesh3d 2 2 2)
add_test(hex-h1-neumann-4x   ${BIN} ${MESHES_DIR}/mesh3d/hex4x.mesh3d 2 2 2)
add_test(hex-h1-neumann-4y   ${BIN} ${MESHES_DIR}/mesh3d/hex4y.mesh3d 2 2 2)
add_test(hex-h1-neumann-4z   ${BIN} ${MESHES_DIR}/mesh3d/hex4z.mesh3d 2 2 2)
add_test(hex-h1-neumann-8    ${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d  2 2 2)

if(FAST_TEST)
	# nothing is here
else(FAST_TEST)
	foreach(a RANGE 2 10)
		foreach(b RANGE 2 10)
			foreach(c RANGE 2 10)
				add_test(hex-h1-neumann-${a}-${b}-${c} ${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d ${a} ${b} ${c})
				add_test(hex-h1-neumann-sc-${a}-${b}-${c} ${BIN} ${MESHES_DIR}/mesh3d/hex1-sc.mesh3d ${a} ${b} ${c})
			endforeach(c)
		endforeach(b)
	endforeach(a)
endif(FAST_TEST)
