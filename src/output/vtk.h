#ifndef _VTK_OUTPUT_ENGINE_H_
#define _VTK_OUTPUT_ENGINE_H_

#include "../output.h"
#include <common/array.h>

/// VTK output engine. 
///
/// TODO: binary file format
///
/// @ingroup visualization
class VtkOutputEngine : public OutputEngine {
public:
	VtkOutputEngine(FILE *file);
	virtual ~VtkOutputEngine();

	/// Run the output with specified output engine
	///
	/// @return true if ok
	/// @param[in] fn A function that will be visualized
	virtual void out(MeshFunction *fn, const char *name, int item = FN_VAL_0);
	virtual void out(Mesh *mesh);

protected:
	void dump_points(MeshFunction *fn);
	void dump_scalars(const char *name, MeshFunction *fn, int item);
	void dump_vectors(const char *name, MeshFunction *fn, int item);
	
	/// file into which the output is done
	FILE *out_file;
	bool has_points;
};

#endif
