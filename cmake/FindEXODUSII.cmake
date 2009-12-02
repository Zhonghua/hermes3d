#
# Exodus2
#
# Looks for library to process exodusII files
# Needs netcdf library installed (this is checked by this module)
#

SET(NETCDF_LIB_SEARCH_PATH
	${NETCDF_ROOT}/lib
	/usr/lib64
	/usr/lib
	/usr/local/lib/
)

SET(EXODUSII_INCLUDE_SEARCH_PATH
	${EXODUSII_ROOT}/include
	/usr/include
	/usr/local/include/
)

SET(EXODUSII_LIB_SEARCH_PATH
	${EXODUSII_ROOT}/lib
	/usr/lib64
	/usr/lib
	/usr/local/lib/
)

FIND_PATH(EXODUSII_INCLUDE_PATH    exodusII.h         ${EXODUSII_INCLUDE_SEARCH_PATH})

FIND_LIBRARY(NETCDF_LIBRARY        netcdf             ${NETCDF_LIB_SEARCH_PATH})
FIND_LIBRARY(EXODUSII_LIBRARY      exoIIv2c           ${EXODUSII_LIB_SEARCH_PATH})


IF(EXODUSII_INCLUDE_PATH)
	SET(EXODUSII_INCLUDE_DIR ${EXODUSII_INCLUDE_DIR} ${EXODUSII_INCLUDE_PATH})
ENDIF(EXODUSII_INCLUDE_PATH)

IF(EXODUSII_LIBRARY)
	SET(EXODUSII_LIBRARIES ${EXODUSII_LIBRARIES} ${EXODUSII_LIBRARY})
	SET(EXODUSII_LIB ${EXODUSII_LIBRARY})
ENDIF(EXODUSII_LIBRARY)

IF(NETCDF_LIBRARY)
	SET(EXODUSII_LIBRARIES ${EXODUSII_LIBRARIES} ${NETCDF_LIBRARY})
ENDIF(NETCDF_LIBRARY)

IF(EXODUSII_INCLUDE_DIR AND EXODUSII_LIBRARIES AND NETCDF_LIBRARY)
	SET(EXODUSII_FOUND TRUE)
ENDIF(EXODUSII_INCLUDE_DIR AND EXODUSII_LIBRARIES AND NETCDF_LIBRARY)

INCLUDE(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EXODUSII DEFAULT_MSG EXODUSII_LIB EXODUSII_INCLUDE_DIR NETCDF_LIBRARY)

