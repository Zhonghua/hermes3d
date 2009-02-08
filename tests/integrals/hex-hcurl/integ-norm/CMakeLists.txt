project(integ-norm)

include(CMake.vars)
include(../../../CMake.common)

add_executable(${PROJECT_NAME}
	main.cc
	${HERMES_COMMON_DIR}/trace.cc
	${HERMES_COMMON_DIR}/timer.cc
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 

include_directories(${hermes3d_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} ${HERMES_CPLX_BIN})

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(integ-norm-1 ${BIN} ${MESHES_DIR}/mesh3d/cube1.mesh3d)
add_test(integ-norm-2 ${BIN} ${MESHES_DIR}/mesh3d/hex1-rot-x90.mesh3d)
add_test(integ-norm-3 ${BIN} ${MESHES_DIR}/mesh3d/hex1-rot-y90.mesh3d)
add_test(integ-norm-4 ${BIN} ${MESHES_DIR}/mesh3d/hex1-rot-z90.mesh3d)

