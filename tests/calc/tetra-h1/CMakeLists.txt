project(calc-tetra-h1)

include(CMake.vars)
include(../../CMake.common)

add_executable(${PROJECT_NAME}
	tetra.cc
	${HERMES_COMMON_DIR}/timer.cc
	${HERMES_COMMON_DIR}/error.cc
	${HERMES_COMMON_DIR}/trace.cc
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
add_test(${PROJECT_NAME}-1  ${BIN} ${MESHES_DIR}/mesh3d/tetra1-ref.mesh3d 4)
#add_test(${PROJECT_NAME}-2  ${BIN} ${MESHES_DIR}/mesh3d/tetra2.mesh3d 3)
#add_test(${PROJECT_NAME}-4  ${BIN} ${MESHES_DIR}/mesh3d/tetra4.mesh3d 2)
#add_test(${PROJECT_NAME}-8  ${BIN} ${MESHES_DIR}/mesh3d/tetra8.mesh3d 1)
