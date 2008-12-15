#ifndef _HDF5_READER_H_
#define _HDF5_READER_H_

#include "meshloader.h"

/// Mesh loader from HDF5 format
///
/// @ingroup meshloaders
class HDF5Reader : public MeshLoader {
public:
	HDF5Reader();
	virtual ~HDF5Reader();

	virtual bool load(const char *file_name, Mesh *mesh);
	virtual bool save(const char *file_name, Mesh *mesh);

	// TODO: save error code and make it accessible via function

	// Mesh attributes
	char *description;					// description of the mesh
};

#endif
