#
# CMake for test the calculations
#

set(FAST_TEST YES)
enable_testing()

if(COMPLEX)
	add_subdirectory(integ-norm)
endif(COMPLEX)
