project(calc-hex-hcurl)

include(CMake.vars)
include(../../CMake.common)

add_executable(${PROJECT_NAME}
	hex-hcurl.cc
	${HERMES_COMMON_DIR}/timer.cc
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

set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(${PROJECT_NAME}-1 ${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d 4)
#add_test(${PROJECT_NAME}-2 ${BIN} ${MESHES_DIR}/mesh3d/cube8-ori.mesh3d 3)
add_test(${PROJECT_NAME}-3 ${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d 3)

