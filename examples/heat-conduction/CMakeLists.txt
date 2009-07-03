project(heat-conduction)

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/CMake.vars)
	include(CMake.vars)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/CMake.vars)


add_executable(${PROJECT_NAME}
	main.cc
	${HERMES_COMMON_DIR}/timer.cc
	${HERMES_COMMON_DIR}/error.cc
)

# Judy
include_directories(${JUDY_INCLUDE_DIR})
include_directories(${hermes3d_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} ${HERMES_REAL_BIN})
target_link_libraries(${PROJECT_NAME} ${JUDY_LIBRARY}) 

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)
