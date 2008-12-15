#
# GLUT
#

FIND_LIBRARY(GLUT_LIBRARY glut ${GLUT_ROOT} /usr/lib64 /usr/lib /usr/local/lib64 /usr/local/lib) 

IF (GLUT_LIBRARY)
	SET(GLUT_FOUND TRUE)
ENDIF (GLUT_LIBRARY)


IF (GLUT_FOUND)
	IF (NOT GLUT_FIND_QUIETLY)
		MESSAGE(STATUS "Found GLUT: ${GLUT_LIBRARY}")
	ENDIF (NOT GLUT_FIND_QUIETLY)
ELSE (GLUT_FOUND)
	IF (GLUT_FIND_REQUIRED)
		MESSAGE(FATAL_ERROR "Could not find GLUT")
	ENDIF (GLUT_FIND_REQUIRED)
ENDIF (GLUT_FOUND)

