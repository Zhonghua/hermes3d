PROJECT_NAME      = Hermes3D
# input
INPUT             = \
	${PROJECT_SOURCE_DIR}/src \
	${PROJECT_SOURCE_DIR}/common
FILE_PATTERNS     = *.cc *.h
RECURSIVE         = YES

# output
OUTPUT_DIRECTORY  = ${PROJECT_BINARY_DIR}/doc
JAVADOC_AUTOBRIEF = YES
QUIET             = YES
WARNINGS          = NO
