#ifndef _GMSH_OUTPUT_ENGINE_H_
#define _GMSH_OUTPUT_ENGINE_H_

#include "../output.h"

/// GMSH output engine.
///
///
///
/// @ingroup visualization
class GmshOutputEngine : public OutputEngine {
public:
	GmshOutputEngine(FILE *file);
	virtual ~GmshOutputEngine();

	/// Run the output with specified output engine
	///
	/// @return true if ok
	/// @param[in] fn A function that will be visualized
	virtual void out(MeshFunction *fn, const char *name, int item = FN_VAL_0);
	virtual void out(Mesh *mesh);

	virtual void out_orders(Space *space, const char *name = "orders");

protected:
	/// file into which the output is done
	FILE *out_file;

	void dump_scalars(int mode, Point3D *pts, double *values, int num_pts);
	void dump_mesh(Mesh *mesh);
};

#endif
