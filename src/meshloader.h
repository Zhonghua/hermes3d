#ifndef _MESHLOADER_H_
#define _MESHLOADER_H_

#include "mesh.h"

/// @defgroup meshloaders Mesh loaders  

/// Abstract class for mesh loaders
///
/// @ingroup meshloaders
class MeshLoader {
public:
	/// Loads the mesh from a file. Aborts the program on error.
	/// @param filename [in] The name of the file.
	/// @param mesh [out] The mesh.
	virtual bool load(const char *file_name, Mesh *mesg) = 0;
	
	/// Saves the mesh, including all refinements, to a file.
	/// Caution: never overwrite hand-created meshes with this function --
	/// all comments in the original file will be lost.
	/// @param filename [in] The name of the file.
	/// @param mesh [out] The mesh.
	virtual bool save(const char *file_name, Mesh *mesh) = 0;
};

#endif
