#
# CMake for tests
#

if(REAL)
	# test cases with REAL version od Hermes3D
	add_subdirectory(mesh3d)
	
	if(USE_HDF5)
	add_subdirectory(hdf5)
	endif(USE_HDF5)
endif(REAL)

