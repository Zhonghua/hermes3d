#
# CMake for examples
#

# HEX
if(WITH_HEX)
    # H1
	if(REAL)
		add_subdirectory(heat-conduction)
	endif(REAL)
endif(WITH_HEX)
