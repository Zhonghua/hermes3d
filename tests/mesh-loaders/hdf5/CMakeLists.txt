project(hdf5-dump)

if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
	include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

add_executable(${PROJECT_NAME}
	hdf5dump.cc
	${HERMES_COMMON_DIR}/trace.cc
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${TESTS_LIBRARIES}) 

include_directories(${hermes3d_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})

CONFIGURE_FILE(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

# Tests

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
# cubes
if(WITH_HEX)
add_test(hdf5-loaders-1 sh -c "${BIN} ${MESHES_DIR}/hdf5/hex1.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/hex1.dump")
add_test(hdf5-loaders-2 sh -c "${BIN} ${MESHES_DIR}/hdf5/hex2.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/hex2.dump")
add_test(hdf5-loaders-3 sh -c "${BIN} ${MESHES_DIR}/hdf5/hex4.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/hex4.dump")
add_test(hdf5-loaders-4 sh -c "${BIN} ${MESHES_DIR}/hdf5/hex8.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/hex8.dump")
endif(WITH_HEX)

# tetras
if(WITH_TETRA)
add_test(hdf5-loaders-tetra-1 sh -c "${BIN} ${MESHES_DIR}/hdf5/tetra1.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/tetra1.dump")
add_test(hdf5-loaders-tetra-2 sh -c "${BIN} ${MESHES_DIR}/hdf5/tetra2.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/tetra2.dump")
add_test(hdf5-loaders-tetra-3 sh -c "${BIN} ${MESHES_DIR}/hdf5/tetra4.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/tetra4.dump")
add_test(hdf5-loaders-tetra-4 sh -c "${BIN} ${MESHES_DIR}/hdf5/tetra8.h5 | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/tetra8.dump")
endif(WITH_TETRA)
