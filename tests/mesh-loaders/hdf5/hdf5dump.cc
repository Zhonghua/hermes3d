/*
 * hdf5dump.cc
 *
 * Tests of mesh loaders
 *
 */

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1
#define ERROR_NOT_ENOUGH_PARAMS						-2

bool testPrint(bool value, const char *msg, bool correct) {
	printf("%s...", msg);
	if (value == correct) {
		printf("OK\n");
		return true;
	}
	else {
		printf("failed\n");
		return false;
	}
}

//
// tests themselves
//

int test_hdf5_loader(char *file_name) {
	Mesh mesh;
	HDF5Reader mloader;
	if (mloader.load(file_name, &mesh)) {
		mesh.dump();
		return ERROR_SUCCESS;
	}
	else {
		printf("failed\n");
		return ERROR_FAILURE;
	}
}

int main(int argc, char *argv[]) {
//	TRACE_START("trace.txt");
	DEBUG_OUTPUT_OFF;
	SET_VERBOSE_LEVEL(0);

	int ret = ERROR_SUCCESS;

	if (argc < 2)
		return ERROR_NOT_ENOUGH_PARAMS;

	if ((ret = test_hdf5_loader(argv[1])) != ERROR_SUCCESS)
		return ret;

	return ret;
}
