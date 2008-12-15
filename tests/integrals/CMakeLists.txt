#
# CMake for testing the integrals
#

enable_testing()

if(REAL)
	if(WITH_HEX)
        add_subdirectory(h1-hex-real)
	endif(WITH_HEX)
endif(REAL)
        