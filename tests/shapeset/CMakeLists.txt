#
# Testing shapesets
#

if(REAL)
	# test cases with REAL version od Hermes3D
	add_subdirectory(hex-h1-lobatto)
	#add_subdirectory(hex-h1-lobatto-ced)
	add_subdirectory(hex-hcurl-lobatto)
endif(REAL)
