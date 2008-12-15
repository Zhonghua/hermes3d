#
# PARDISO
#

FIND_LIBRARY(PARDISO_LIBRARY ${PARDISO_LIB} ${PARDISO_ROOT}) 

IF (PARDISO_LIBRARY)
	SET(PARDISO_FOUND TRUE)
ENDIF (PARDISO_LIBRARY)


IF (PARDISO_FOUND)
	IF (NOT PARDISO_FIND_QUIETLY)
		MESSAGE(STATUS "Found pardiso: ${PARDISO_LIBRARY}")
	ENDIF (NOT PARDISO_FIND_QUIETLY)
ELSE (PARDISO_FOUND)
	IF (PARDISO_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could not find pardiso")
	ENDIF (PARDISO_FIND_REQUIRED)
ENDIF (PARDISO_FOUND)

