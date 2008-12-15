/*
 * electordes.cc
 *
 * - \laplace u = 0
 * na hranici c. 1, 7, 8, 9 je Dirichlet, u = 0
 * na hranici c. 2, 10, 11, 12 je Dirichlet, u = 10000
 * na hranici c. 3, 4, 5, 6 je Neumann, du/dv = 0
 *
 * rad elementu: 1, prip. 2
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

double fnc(double x, double y, double z) {
	return 0.0;
}

double dfnc(double x, double y, double z) {
	return 0.0;
}

double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 0.0;
	dy = 0.0;
	dz = 0.0;

	return 0.0;
}

////

EBCType bc_types(int marker) {
	switch (marker) {
		case 1:
		case 7:
		case 8:
		case 9:
		case 2:
		case 10:
		case 11:
		case 12:
			return BC_ESSENTIAL;
		case 3:
		case 4:
		case 5:
		case 6:
			return BC_NATURAL;
		default:
			EXIT(ERR_FAILURE, "Invalid marker.");
	}
	return BC_NATURAL;
}

double bc_values(int marker, double x, double y, double z) {
	switch (marker) {
		case 1:
		case 7:
		case 8:
		case 9:
			return 0;
		case 2:
		case 10:
		case 11:
		case 12:
			return 10000;
		case 3:
		case 4:
		case 5:
		case 6:
			return 0;
		default: EXIT(ERR_FAILURE, "Unknown marker.");
	}
}

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	return int_F_v(dfnc, fv, rv);
}

scalar linear_form_surf(RealFunction *fv, RefMap *rv, FacePos *fp) {
	return surf_int_G_v(fv, rv, fp);
}

// output ///
void dump_output(Solution &sln) {
#ifdef OUTPUT_DIR
//	if (argv[3] == NULL || strcasecmp(argv[3], "dump") != 0) {
		printf("* Dumping output...\n");
		// output
		char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			GmshOutputEngine output(ofile);
			output.out(&sln, "Uh", FN_VAL_0);
//			output.out(&sln, "Uh dx", FN_DX_0, 0);
//			output.out(&sln, "Uh dy", FN_DY_0, 0);
//			output.out(&sln, "Uh dz", FN_DZ_0, 0);

			fclose(ofile);
		}
		else {
			ERROR("Can't open '%s' for writing.", of_name);
		}
//	}
#endif
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
	Mesh3DReader mesh_loader;
	Mesh mesh;
	if (!mesh_loader.load(argv[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", argv[1]);
		return ERR_FAILURE;
	}

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_bc_values(bc_values);

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
	d.set_linear_form(0, linear_form, linear_form_surf);

	// assemble siffness matrix
	d.create_stiffness_matrix();
	printf("  - total matrix size: %0.1lf MB\n", (double) solver.get_matrix_size() / (1024*1024));

	Timer assemble_timer("assembling stiffness matrix");
	assemble_timer.start();
	d.assemble_stiffness_matrix_and_rhs();
	assemble_timer.stop();
	printf("  - %s: %s (%lf secs)\n", assemble_timer.get_name(), assemble_timer.get_human_time(), assemble_timer.get_seconds());

	Solution sln(&mesh);

	if (argv[3] != NULL) {
		if (strcasecmp(argv[3], "dump") == 0) {
			// dump <data file>
	//		const char *m_file_name = "a.m";
	//		FILE *mfile = fopen(m_file_name, "w");
			FILE *mfile = fopen(argv[4], "w");
			if (mfile != NULL) {
				solver.dump_matrix(mfile, "A", DF_HERMES_BIN);
				solver.dump_rhs(mfile, "b", DF_HERMES_BIN);
				
				fclose(mfile);
			}
			else {
				ERROR("Can not open file '%s' for writing", argv[4]);
			}
		}
		else if (strcasecmp(argv[3], "load") == 0) {
			// load <data file>

//			Solution sln(mesh);
			int n = ndofs;	// 313660;
//			swprintf(argv[5], "%d", &n);
			double *s = new double[n + 1];
			s[0] = 1.0;

			FILE *sfile = fopen(argv[4], "r");
			for (int i = 1; i <= n; i++) {
				char str[256];
				fgets(str, 255, sfile);
				double d;
				sscanf(str, "%lf", &d);
				s[i] = d;
			}
			fclose(sfile);

			sln.set_space_and_pss(&space, &pss);
			sln.set_solution_vector(s, false);
			
			dump_output(sln);
			
			delete [] s;
		}
	}
	else {
		// solve the stiffness matrix
		Timer solve_timer("Solving stiffness matrix");
		solve_timer.start();
		bool solved = d.solve_system(1, &sln);
		solve_timer.stop();
		// output the measured values
		printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(), solve_timer.get_seconds());
		
		dump_output(sln);
	}

/*	FILE *sfile = fopen("sln", "w");
	double *s = sln.get_solution_vector();
	for (int i = 0; i <= ndofs; i++) {
		fprintf(sfile, "%.18e\n", s[i]);
	}
*/

//	if (solved) {
/*		printf("* Solution:\n");
		double *s = sln.get_solution_vector();
		for (int i = 1; i <= ndofs; i++) {
			printf(" x[% 3d] = % lf\n", i, s[i]);
		}
*/

//	}


#ifdef USE_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

