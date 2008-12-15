project(hex-hcurl)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	hex-hcurl.cc
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
add_test(calc-hex-hcurl-1 ${BIN} ${MESHES_DIR}/mesh3d/cube1.mesh3d 4)
add_test(calc-hex-hcurl-2 ${BIN} ${MESHES_DIR}/mesh3d/cube8-ori.mesh3d 3)
add_test(calc-hex-hcurl-3 ${BIN} ${MESHES_DIR}/mesh3d/cube8.mesh3d 3)

