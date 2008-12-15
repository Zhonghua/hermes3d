project(mesh3d-dump)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	mesh3ddump.cc
	${HERMES_COMMON_DIR}/trace.cc
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

set(BIN ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})

# cubes
if(WITH_HEX)
add_test(m3d-loaders-1 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/hex1.dump")
add_test(m3d-loaders-2 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/hex2.dump")
add_test(m3d-loaders-3 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/hex4.dump")
add_test(m3d-loaders-4 sh -c "${BIN} ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/hex8.dump")

add_test(m3d-loaders-5 sh -c "${BIN} ${MESHES_DIR}/mesh3d/missing_boundary.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/missing_boundary.dump")
add_test(m3d-loaders-6 sh -c "${BIN} ${MESHES_DIR}/mesh3d/invalid_vertex_indices.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/invalid_vertex_indices.dump")
add_test(m3d-loaders-7 sh -c "${BIN} ${MESHES_DIR}/mesh3d/missing_outer_face_boundary.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/missing_outer_face_boundary.dump")

endif(WITH_HEX)

# tetras
if(WITH_TETRA)
add_test(m3d-loaders-tetra-1 sh -c "${BIN} ${MESHES_DIR}/mesh3d/tetra1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/tetra1.dump")
add_test(m3d-loaders-tetra-2 sh -c "${BIN} ${MESHES_DIR}/mesh3d/tetra2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/tetra2.dump")
add_test(m3d-loaders-tetra-3 sh -c "${BIN} ${MESHES_DIR}/mesh3d/tetra4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/tetra4.dump")
add_test(m3d-loaders-tetra-4 sh -c "${BIN} ${MESHES_DIR}/mesh3d/tetra8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/../dump/tetra8.dump")
endif(WITH_TETRA)


