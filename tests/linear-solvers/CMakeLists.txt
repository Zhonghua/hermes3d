project(solver)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	main.cc
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

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
if(USE_PETSC)
	add_test(petsc-solver-1 sh -c "${BIN} petsc ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(petsc-solver-2 sh -c "${BIN} petsc ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
endif(USE_PETSC)

if(USE_UMFPACK)
	add_test(umfpack-solver-1 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(umfpack-solver-2 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(umfpack-solver-3 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	add_test(umfpack-solver-4 sh -c "${BIN} umfpack ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
	
	add_test(umfpack-solver-b-1 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(umfpack-solver-b-2 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(umfpack-solver-b-3 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	add_test(umfpack-solver-b-4 sh -c "${BIN} umfpack-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
endif(USE_UMFPACK)

if(USE_PARDISO)
	add_test(pardiso-solver-1 sh -c "${BIN} pardiso ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(pardiso-solver-2 sh -c "${BIN} pardiso ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(pardiso-solver-3 sh -c "${BIN} pardiso ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	add_test(pardiso-solver-4 sh -c "${BIN} pardiso ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
	
	add_test(pardiso-solver-b-1 sh -c "${BIN} pardiso-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-1 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-1")
	add_test(pardiso-solver-b-2 sh -c "${BIN} pardiso-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-2 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-2")
	add_test(pardiso-solver-b-3 sh -c "${BIN} pardiso-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-3 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-3")
	add_test(pardiso-solver-b-4 sh -c "${BIN} pardiso-block ${CMAKE_CURRENT_SOURCE_DIR}/in/linsys-singular | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/linsys-singular")
endif(USE_PARDISO)

