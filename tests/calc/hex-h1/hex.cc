/*
 * hex.cc
 *
 *
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

// first two Lobatto shape functions
#define l0(x) ((1.0 - (x)) * 0.5)
#define l1(x) ((1.0 + (x)) * 0.5)

// error should be smaller than this epsilon
#define EPS								10e-10F

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = -0.5 * x * l0(y) * l1(y) * l0(z) * l1(z);
	dy = -0.5 * y * l0(x) * l1(x) * l0(z) * l1(z);
	dz = -0.5 * z * l0(x) * l1(x) * l0(y) * l1(y);

	return l0(x) * l1(x) * l0(y) * l1(y) * l0(z) * l1(z);
}

//

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

double f(double x, double y, double z) {
	return
		0.5 * l0(y) * l1(y) * l0(z) * l1(z) +
		0.5 * l0(x) * l1(x) * l0(z) * l1(z) +
		0.5 * l0(x) * l1(x) * l0(y) * l1(y);
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	return int_F_v(f, fv, rv);
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

	if (argc < 3) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	H1ShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[1]);
		return ERR_FAILURE;
	}

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);

	int o;
	sscanf(args[2], "%d", &o);
	order3_t order(o, o, o);
	printf("  - Setting uniform order to %s\n", order.str());
	space.set_uniform_order(order);

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
	Solution sln(&mesh);
	bool solved = d.solve_system(1, &sln);
	solve_timer.stop();

	// output the measured values
	printf("%s: %s (%lf secs)\n", assemble_timer.get_name(), assemble_timer.get_human_time(), assemble_timer.get_seconds());
	printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(), solve_timer.get_seconds());

	if (solved) {
		printf("* Solution:\n");
		double *s = sln.get_solution_vector();
		for (int i = 0; i <= ndofs; i++) {
			printf(" x[% 3d] = % lf\n", i, s[i]);
		}

//		s[1] = 4.166667e-02;
//		s[2] = s[3] = -1.701035e-02;

		ExactSolution ex_sln(&mesh, exact_solution);

		// norm
		double h1_sln_norm = h1_norm(&sln);
		double h1_err_norm = h1_error(&sln, &ex_sln);
		printf(" - H1 solution norm:   % le\n", h1_sln_norm);
		printf(" - H1 error norm:      % le\n", h1_err_norm);

		double l2_sln_norm = l2_norm(&sln);
		double l2_err_norm = l2_error(&sln, &ex_sln);
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

			GmshOutputEngine out_eng(ofile);
			Output output(&out_eng);
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
			ERROR("Can not not open '%s' for writing.", of_name);
		}
#endif
	}

#ifdef USE_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

