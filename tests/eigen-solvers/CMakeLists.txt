project(eigen-solver)

if(HAVE_EIGEN_SOLVER)

if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
    include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

add_executable(${PROJECT_NAME}
	main.cc
	${HERMES_COMMON_DIR}/trace.cc
)

target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 
target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
if(WITH_SLEPC)
	add_test(eigen-solver-slepc-1 sh -c "${BIN} slepc | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/slepc")
	add_test(eigen-solver-slepc-2 sh -c "${BIN} slepc-gen | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/slepc-gen")
endif(WITH_SLEPC)

endif(HAVE_EIGEN_SOLVER)
