#ifndef _MESH3DLOADER_H_
#define _MESH3DLOADER_H_

#include "../meshloader.h"

/// Mesh loader from Mesh3D format
///
/// @ingroup meshloaders
class Mesh3DReader : public MeshLoader {
public:
	Mesh3DReader();
	virtual ~Mesh3DReader();

	virtual bool load(const char *file_name, Mesh *mesh);
	virtual bool save(const char *file_name, Mesh *mesh);
};

#endif
