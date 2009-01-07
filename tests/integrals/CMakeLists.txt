#
# CMake for testing the integrals
#

enable_testing()

if(REAL)
	if(WITH_HEX)
        add_subdirectory(hex-h1-real)
	endif(WITH_HEX)
endif(REAL)
        