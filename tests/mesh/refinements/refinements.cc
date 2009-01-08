/*
 * refinements.cc
 *
 * usage: $0 <mesh file> <element id> <refinement id> [<element id> <refinement id>...]
 *
 */

#include "config.h"
#ifdef USE_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>


// helpers ////////////////////////////////////////////////////////////////////////////////////////

int parse_reft(char *str) {
	if (strcasecmp(str, "x") == 0) return REFT_HEX_X;
	else if (strcasecmp(str, "y") == 0) return REFT_HEX_Y;
	else if (strcasecmp(str, "z") == 0) return REFT_HEX_Z;
	else if (strcasecmp(str, "xy") == 0 || strcasecmp(str, "yx") == 0) return REFT_HEX_XY;
	else if (strcasecmp(str, "xz") == 0 || strcasecmp(str, "zx") == 0) return REFT_HEX_XZ;
	else if (strcasecmp(str, "yz") == 0 || strcasecmp(str, "zy") == 0) return REFT_HEX_YZ;
	else if (strcasecmp(str, "xyz") == 0) return REFT_HEX_XYZ;
	else return REFT_HEX_NONE;
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;

#ifdef USE_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	if (argc < 1) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[1]);
		return ERR_FAILURE;
	}

	// apply refinements
	bool ok = true;
	for (int k = 2; k < argc && ok; k += 2) {
		int elem_id, reft_id;
		sscanf(args[k], "%d", &elem_id);
		reft_id = parse_reft(args[k + 1]);
		ok = mesh.refine_element(elem_id, reft_id);
	}

	if (ok) {
		mesh.dump();
	}
	else {
		ERROR("Unable to refine a mesh.");
		res = ERR_FAILURE;
	}


#ifdef USE_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}
