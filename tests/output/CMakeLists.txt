project(output)

if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
	include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

set(BIN_GMSH ${PROJECT_NAME}-gmsh)
set(BIN_VTK  ${PROJECT_NAME}-vtk)


include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${TESTS_INCLUDE_DIRS})
include_directories(${hermes3d_SOURCE_DIR})

# GMSH ####

add_executable(${BIN_GMSH}
	main.cc
)

set_target_properties(${BIN_GMSH}
    PROPERTIES
	COMPILE_FLAGS -DGMSH)


target_link_libraries(${BIN_GMSH} ${TESTS_LIBRARIES}) 
target_link_libraries(${BIN_GMSH} ${HERMES_REAL_BIN})

# Tests

# GMSH ####
set(BIN ${PROJECT_BINARY_DIR}/${BIN_GMSH})

# solution
add_test(${PROJECT_NAME}-gmsh-sln-1 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-1")
add_test(${PROJECT_NAME}-gmsh-sln-2 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-2")
add_test(${PROJECT_NAME}-gmsh-sln-3 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-4")
add_test(${PROJECT_NAME}-gmsh-sln-4 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-8")
add_test(${PROJECT_NAME}-gmsh-sln-5 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-27")
add_test(${PROJECT_NAME}-gmsh-sln-6 sh -c "${BIN} sln ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-sln-fichera-corner")

# order
add_test(${PROJECT_NAME}-gmsh-ord-1 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex1.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-1")
add_test(${PROJECT_NAME}-gmsh-ord-2 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex2.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-2")
add_test(${PROJECT_NAME}-gmsh-ord-3 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex4.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-4")
add_test(${PROJECT_NAME}-gmsh-ord-4 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex8.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-8")
add_test(${PROJECT_NAME}-gmsh-ord-5 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/hex27.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-27")
add_test(${PROJECT_NAME}-gmsh-ord-6 sh -c "${BIN} ord ${MESHES_DIR}/mesh3d/fichera-corner.mesh3d | diff - ${CMAKE_CURRENT_SOURCE_DIR}/out/gmsh-ord-fichera-corner")

