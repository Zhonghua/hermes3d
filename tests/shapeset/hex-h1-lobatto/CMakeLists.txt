project(hex-h1-lobatto-shapeset)

include(CMake.vars)

if(WITH_HEX)
#
add_executable(${PROJECT_NAME}
	main.cc
	lindep.cc
	zero.cc
	cont.cc
	grad.cc
	graddir.cc
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
add_test(hex-h1-lobbatto-shapeset-1 ${BIN})

endif(WITH_HEX)

