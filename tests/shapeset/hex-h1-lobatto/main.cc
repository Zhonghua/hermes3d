/*
 * main.cc
 *
 * Test for H1 lobatto shapeset for Hex
 *
 */

#include "config.h"
#ifdef USE_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

// forward declarations (simpler than defining a special header file for each module that exports just one function)
bool test_lin_indep(Shapeset *shapeset);
bool test_zero_values(Shapeset *shapeset);
bool test_continuity(Shapeset *shapeset);
bool test_gradients(Shapeset *shapeset);
bool test_gradients_directly(Shapeset *shapeset);

//
// main
//
int main(int argc, char *argv[]) {
	int res = ERR_SUCCESS;

#ifdef USE_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif

	H1ShapesetLobattoHex shapeset;

	try {
		// I. linear independency
		if (!test_lin_indep(&shapeset))
			throw ERR_FAILURE;

		// II. test zero fn. values
		if (!test_zero_values(&shapeset))
			throw ERR_FAILURE;

		// III. continuity on boundaries
		if (!test_continuity(&shapeset))
			throw ERR_FAILURE;

		// IV. gradients
		if (!test_gradients(&shapeset))
			throw ERR_FAILURE;

		// V. computes gradients numericaly from fn values and compares
		if (!test_gradients_directly(&shapeset))
			throw ERR_FAILURE;

		printf("Shapeset OK\n");
	}
	catch (int e) {
		printf("Failed\n");
		res = e;
	}

#ifdef USE_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}
