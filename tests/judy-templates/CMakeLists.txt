project(judy-templates)

if(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)
	include(CMake.vars)
endif(EXISTS ${PROJECT_SOURCE_DIR}/CMake.vars)

add_executable(${PROJECT_NAME}
	main.cc
	${HERMES_COMMON_DIR}/trace.cc
)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Judy
include_directories(${JUDY_INCLUDE_DIR})
target_link_libraries(${PROJECT_NAME} ${JUDY_LIBRARY}) 

# Tests

add_test(judy-templates-1 ${PROJECT_NAME})

