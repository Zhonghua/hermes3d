project(timer)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	main.cc
	${HERMES_COMMON_DIR}/trace.cc
	${HERMES_COMMON_DIR}/timer.cc
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 

include_directories(${hermes3d_SOURCE_DIR})
if(REAL)
	target_link_libraries(${PROJECT_NAME} ${HERMES_BIN})
elseif(COMPLEX)
	target_link_libraries(${PROJECT_NAME} ${HERMES_CPLX_BIN})
endif(REAL)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

#
# Tests
#

enable_testing()

add_test(timer-1 ${PROJECT_NAME})

