project(qorder-solution)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	main.cc
	${HERMES_COMMON_DIR}/trace.cc
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

add_test(qorder-solution-1 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex1.mesh3d)
add_test(qorder-solution-2 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex1-sc.mesh3d)
add_test(qorder-solution-3 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex2.mesh3d)
add_test(qorder-solution-4 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex4.mesh3d)
add_test(qorder-solution-5 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex8.mesh3d)
add_test(qorder-solution-6 ${PROJECT_NAME} ${MESHES_DIR}/mesh3d/hex27.mesh3d)
