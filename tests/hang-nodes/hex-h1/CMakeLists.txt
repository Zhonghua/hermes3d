project(hnnd-hex-h1)

include(CMake.vars)

set(BIN_DIRICHLET ${PROJECT_NAME}-dirichlet)
set(BIN_NEWTON ${PROJECT_NAME}-newton)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIR})
include_directories(${hermes3d_SOURCE_DIR})

add_executable(${BIN_DIRICHLET}
	hang-nodes-cont.cc
	${HERMES_COMMON_DIR}/timer.cc
	${HERMES_COMMON_DIR}/error.cc
)

target_link_libraries(${BIN_DIRICHLET} ${TESTS_LIBRARIES}) 
target_link_libraries(${BIN_DIRICHLET} ${HERMES_BIN})
set_target_properties(${BIN_DIRICHLET}
    PROPERTIES
	COMPILE_FLAGS -DDIRICHLET)

add_executable(${BIN_NEWTON}
	hang-nodes-cont.cc
	${HERMES_COMMON_DIR}/timer.cc
	${HERMES_COMMON_DIR}/error.cc
)
target_link_libraries(${BIN_NEWTON} ${TESTS_LIBRARIES}) 
target_link_libraries(${BIN_NEWTON} ${HERMES_BIN})
set_target_properties(${BIN_NEWTON}
    PROPERTIES
	COMPILE_FLAGS -DNEWTON)


configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

#
# Tests
#

enable_testing()

set(BIN ${PROJECT_BINARY_DIR}/${BIN_DIRICHLET})
set(NAME ${PROJECT_NAME}-dir)
include(../CMake.tests.hex)

set(BIN ${PROJECT_BINARY_DIR}/${BIN_NEWTON})
set(NAME ${PROJECT_NAME}-nwt)
include(../CMake.tests.hex)

