project(judy-templates)

include(CMake.vars)

add_executable(${PROJECT_NAME}
	main.cc
	${HERMES_COMMON_DIR}/trace.cc
)

configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake
	${CMAKE_CURRENT_SOURCE_DIR}/config.h
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_package(JUDY REQUIRED)
if(JUDY_FOUND)
	include_directories(${JUDY_INCLUDE_DIR})
	target_link_libraries(${PROJECT_NAME} ${JUDY_LIBRARY}) 
endif(JUDY_FOUND)

#
# Tests
#

enable_testing()

add_test(judy-templates-1 ${PROJECT_NAME})

