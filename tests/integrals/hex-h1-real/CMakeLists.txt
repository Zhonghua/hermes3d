project(intgrl-hex-h1-real)

include(CMake.vars)
include(../../CMake.common)

add_executable(${PROJECT_NAME}
	h1-hex-real.cc
	${HERMES_COMMON_DIR}/trace.cc
	${HERMES_COMMON_DIR}/error.cc
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 

include_directories(${hermes3d_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

add_test(${PROJECT_NAME}-1 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex1.mesh3d)
add_test(${PROJECT_NAME}-2 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex2.mesh3d)
add_test(${PROJECT_NAME}-3 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex4.mesh3d)
add_test(${PROJECT_NAME}-4 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex8.mesh3d)

