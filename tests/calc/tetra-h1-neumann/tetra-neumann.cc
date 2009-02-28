// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

/*
 * tetra-neumann.cc
 *
 *
 *
 */

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

// error should be smaller than this epsilon
#define EPS								10e-10F

#define ERR_LOADING_MESH				-1


#define STARTING_DOF 0
#define OUTPUT_PRECISION 1

double fnc(double x, double y, double z) {
	return x*x + y*y + z*z;
}

double dfnc(double x, double y, double z) {
	return -6.0 + fnc(x, y, z);
}

double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return fnc(x, y, z);
}

////

EBCType bc_types(int marker) {
	return BC_NATURAL;
}

double bc_values(int marker, double x, double y, double z, int comp) {
	switch (marker) {
		case 1: return -2 * y;
		case 2: return (2.0 / 7.0) * (3 * x + 2 * y + 6 * z);
		case 3: return -2 * x;
		case 4: return -2 * z;
		default: EXIT(ERR_FAILURE, "Unknown marker"); return 0.0;
	}
}

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv) + int_u_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	return int_F_v(dfnc, fv, rv);
}

scalar linear_form_surf(RealFunction *fv, RefMap *rv, FacePos *fp) {
	return surf_int_G_v(fv, rv, fp);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
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
	space.set_bc_values(bc_values);

	int o;
	sscanf(argv[2], "%d", &o);
	printf("  - Setting uniform order to %d\n", o);
	space.set_uniform_order(o);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

#if defined WITH_UMFPACK
	UMFPackLinearSolver solver;
#elif defined WITH_PARDISO
	PardisoLinearSolver solver;
#elif defined WITH_PETSC
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

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

