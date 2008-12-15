/*
 * tetra.cc
 *
 *
 *
 */

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>
#ifdef USE_PETSC
#include <petsc.h>
#endif

// error should be smaller than this epsilon
#define EPS								10e-10F

#define ERR_LOADING_MESH				-1


#define STARTING_DOF 0
#define OUTPUT_PRECISION 1

#define la0(x,y,z) (((y) + 1) / 2)
#define la1(x,y,z) (-(1 + (x) + (y) + (z)) / 2)
#define la2(x,y,z) (((x) + 1) / 2)
#define la3(x,y,z) (((z) + 1) / 2)
#define ph0 (-2.0 * 1.22474487139158904909864203735)

double fnc(double x, double y, double z) {
	return la0(x, y, z) * la1(x, y, z) * la2(x, y, z) * la3(x, y, z) * ph0 * ph0 * ph0;
}

double dfnc(double x, double y, double z) {
	return
		- 1.5 * sqrt(3.0 / 2.0) * (y + 1) * (z + 1)
		- 1.5 * sqrt(3.0 / 2.0) * (x + 1) * (z + 1)
		- 1.5 * sqrt(3.0 / 2.0) * (x + 1) * (y + 1);
}

double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((y + 1) * (- z - y - x - 1) * (z + 1)));
	dy = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((x + 1) * (- z - y - x - 1) * (z + 1)));
	dz = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((x + 1) * (y + 1) * (- z - y - x - 1)));

	return fnc(x, y, z);
}

////

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction *fu, RefMap *ru) {
	return int_F_u(dfnc, fu, ru);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	int res = ERR_SUCCESS;

#ifdef USE_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	if (argc < 3) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	H1ShapesetLobattoTetra shapeset;
	PrecalcShapeset pss(&shapeset);

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", argv[1]);
		return ERR_FAILURE;
	}

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);

	int o;
	sscanf(argv[2], "%d", &o);
	printf("  - Setting uniform order to %d\n", o);
	space.set_uniform_order(o);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

#if defined USE_UMFPACK
	UMFPackLinearSolver solver;
#elif defined USE_PARDISO
	PardisoLinearSolver solver;
#elif defined USE_PETSC
	PetscLinearSolver solver;
#endif

	Discretization d(&solver);
	d.set_num_equations(1);
	d.set_spaces(1, &space);
	d.set_pss(1, &pss);

	d.set_bilinear_form(0, 0, bilinear_form);
	d.set_linear_form(0, linear_form);

	// assemble siffness matrix
	d.create_stiffness_matrix();

	Timer assemble_timer("Assembling stiffness matrix");
	assemble_timer.start();
	d.assemble_stiffness_matrix_and_rhs();
	assemble_timer.stop();

	// solve the stiffness matrix
	Timer solve_timer("Solving stiffness matrix");
	solve_timer.start();
	Solution sln(mesh);
	bool solved = d.solve_system(1, &sln);
	solve_timer.stop();

	// output the measured values
	printf("%s: %s (%lf secs)\n", assemble_timer.get_name(), assemble_timer.get_human_time(), assemble_timer.get_seconds());
	printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(), solve_timer.get_seconds());

	if (solved) {
/*		printf("* Solution:\n");
		double *s = sln.get_solution_vector();
		for (int i = 1; i <= ndofs; i++) {
			printf(" x[% 3d] = % lf\n", i, s[i]);
		}
*/
		// norm
		double h1_sln_norm = h1_norm(&sln);
		double h1_err_norm = h1_error_norm_exact(&sln, exact_solution);

		printf(" - H1 solution norm:   % le\n", h1_sln_norm);
		printf(" - H1 error norm:      % le\n", h1_err_norm);

		double l2_sln_norm = l2_norm(&sln);
		double l2_err_norm = l2_error_norm_exact(&sln, exact_solution);
		printf(" - L2 solution norm:   % le\n", l2_sln_norm);
		printf(" - L2 error norm:      % le\n", l2_err_norm);

		if (h1_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

#ifdef OUTPUT_DIR
		// output
		char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			ExactSolution ex_sln(&mesh, exact_solution);
			DiffFilter eh(&mesh, &sln, &ex_sln);
//			DiffFilter eh_dx(&mesh, &sln, &ex_sln, FN_DX, FN_DX);
//			DiffFilter eh_dy(&mesh, &sln, &ex_sln, FN_DY, FN_DY);
//			DiffFilter eh_dz(&mesh, &sln, &ex_sln, FN_DZ, FN_DZ);
						
			GmshOutputEngine output(ofile);
			output.out(&sln, "Uh");
//			output.out(&sln, "Uh dx", FN_DX_0);
//			output.out(&sln, "Uh dy", FN_DY_0);
//			output.out(&sln, "Uh dz", FN_DZ_0);
			output.out(&eh, "Eh");
//			output.out(&eh_dx, "Eh dx");
//			output.out(&eh_dy, "Eh dy");
//			output.out(&eh_dz, "Eh dz");
			output.out(&ex_sln, "U");
//			output.out(&ex_sln, "U dx", FN_DX_0);
//			output.out(&ex_sln, "U dy", FN_DY_0);
//			output.out(&ex_sln, "U dz", FN_DZ_0);

			fclose(ofile);
		}
		else {
			ERROR("Can't open '%s' for writing.", of_name);
		}
#endif
	}

#ifdef USE_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

