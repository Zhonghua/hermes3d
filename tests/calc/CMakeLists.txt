#
# CMake for test the calculations
#

# TETRA
if(WITH_TETRA)
    if(REAL)
    	add_subdirectory(tetra-h1)
    	add_subdirectory(tetra-h1-dirichlet)
    	add_subdirectory(tetra-h1-neumann-ref)
    	add_subdirectory(tetra-h1-neumann)
    	add_subdirectory(tetra-h1-newton-ref)
    	add_subdirectory(tetra-h1-newton)
    endif(REAL)
endif(WITH_TETRA)

# HEX
if(WITH_HEX)
    if(REAL)
	add_subdirectory(hex-h1)
	add_subdirectory(hex-h1-1face)
	add_subdirectory(hex-h1-dirichlet)
	add_subdirectory(hex-h1-neumann)
	add_subdirectory(hex-h1-newton)
    endif(REAL)
endif(WITH_HEX)
