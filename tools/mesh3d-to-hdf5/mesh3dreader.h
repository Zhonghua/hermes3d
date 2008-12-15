#ifndef MESH3D_READER_H_
#define MESH3D_READER_H_

class Mesh;

/// Mesh loader from Mesh3D format
///
/// @ingroup meshloaders
class Mesh3DLoader {
public:
	Mesh3DLoader();
	virtual ~Mesh3DLoader();

	virtual bool load(const char *file_name, Mesh *mesh);
};

#endif /*MESH3D_READER_H_*/
