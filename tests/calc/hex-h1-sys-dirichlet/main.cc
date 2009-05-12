// This file is part of Hermes3D
//
// Copyright (c) 2009 David Andrs <dandrs@unr.edu>
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

// Two uncoupled equations:
// -\laplace u_1 = f_1        u_1 = x^2 + y^2 + z^2)
// -\laplace u_2 = f_2        u_2 = x^3 + y^3 + z^3
//
// u_1 = 0 on \partial\Omega
// u_2 = 0 on \partial\Omega
//

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

double u1(double x, double y, double z) {
	return x*x + y*y + z*z;
}

double u2(double x, double y, double z) {
	return x*x*x + y*y*y + z*z*z;
}

// needed for calculation norms and used by visualizator
double exact_sln_fn_1(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return u1(x, y, z);
}

double exact_sln_fn_2(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 3 * x*x;
	dy = 3 * y*y;
	dz = 3 * z*z;

	return u2(x, y, z);
}

//

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

double bc_values_1(int marker, double x, double y, double z) {
	return u1(x, y, z);
}

double bc_values_2(int marker, double x, double y, double z) {
	return u2(x, y, z);
}

scalar bilinear_form_1(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar bilinear_form_2(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

double f1(double x, double y, double z) {
	return -6.0;
}

scalar linear_form_1(RealFunction *fv, RefMap *rv) {
	return int_F_v(f1, fv, rv);
}

double f2(double x, double y, double z) {
	return -(6 * x + 6 * y + 6 * z);
}

scalar linear_form_2(RealFunction *fv, RefMap *rv) {
	return int_F_v(f2, fv, rv);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	if (argc < 2) {
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

	printf("* Setup space #1\n");
	H1Space space1(&mesh, &shapeset);
	space1.set_bc_types(bc_types);
	space1.set_bc_values(bc_values_1);

	order3_t o1(2, 2, 2);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o1.x, o1.y, o1.z);
	space1.set_uniform_order(o1);

	printf("* Setup space #2\n");
	H1Space space2(&mesh, &shapeset);
	space2.set_bc_types(bc_types);
	space2.set_bc_values(bc_values_2);

	order3_t o2(3, 3, 3);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o2.x, o2.y, o2.z);
	space2.set_uniform_order(o2);

	int ndofs = 0;
	ndofs += space1.assign_dofs();
	ndofs += space2.assign_dofs(ndofs);
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

#if defined WITH_UMFPACK
	UMFPackMatrix mat;
	UMFPackVector rhs;
	UMFPackLinearSolver solver(mat, rhs);
#elif defined WITH_PARDISO
	PardisoLinearSolver solver;
#elif defined WITH_PETSC
	PetscMatrix mat;
	PetscVector rhs;
	PetscLinearSolver solver(mat, rhs);
#endif

	Discretization d;
	d.set_num_equations(2);
	d.set_spaces(2, &space1, &space2);
	d.set_pss(1, &pss);

	d.set_bilinear_form(0, 0, bilinear_form_1);
	d.set_linear_form(0, linear_form_1);
	d.set_bilinear_form(1, 1, bilinear_form_2);
	d.set_linear_form(1, linear_form_2);

	// assemble stiffness matrix
	d.create(&mat, &rhs);

	Timer assemble_timer("Assembling stiffness matrix");
	assemble_timer.start();
	d.assemble(&mat, &rhs);
	assemble_timer.stop();

	// solve the stiffness matrix
	Timer solve_timer("Solving stiffness matrix");
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();

	// output the measured values
	printf("%s: %s (%lf secs)\n", assemble_timer.get_name(), assemble_timer.get_human_time(), assemble_timer.get_seconds());
	printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(), solve_timer.get_seconds());

	if (solved) {
		// solution 1
		Solution sln1(&mesh);
		sln1.set_space_and_pss(&space1, &pss);
		sln1.set_solution_vector(solver.get_solution(), false);

		ExactSolution esln1(&mesh, exact_sln_fn_1);
		// norm
		double h1_sln_norm1 = h1_norm(&sln1);
		double h1_err_norm1 = h1_error(&sln1, &esln1);

		printf(" - H1 solution norm:   % le\n", h1_sln_norm1);
		printf(" - H1 error norm:      % le\n", h1_err_norm1);

		double l2_sln_norm1 = l2_norm(&sln1);
		double l2_err_norm1 = l2_error(&sln1, &esln1);
		printf(" - L2 solution norm:   % le\n", l2_sln_norm1);
		printf(" - L2 error norm:      % le\n", l2_err_norm1);

		if (h1_err_norm1 > EPS || l2_err_norm1 > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

		// solution 2
		Solution sln2(&mesh);
		sln2.set_space_and_pss(&space2, &pss);
		sln2.set_solution_vector(solver.get_solution(), false);

		ExactSolution esln2(&mesh, exact_sln_fn_2);
		// norm
		double h1_sln_norm2 = h1_norm(&sln2);
		double h1_err_norm2 = h1_error(&sln2, &esln2);

		printf(" - H1 solution norm:   % le\n", h1_sln_norm2);
		printf(" - H1 error norm:      % le\n", h1_err_norm2);

		double l2_sln_norm2 = l2_norm(&sln2);
		double l2_err_norm2 = l2_error(&sln2, &esln2);
		printf(" - L2 solution norm:   % le\n", l2_sln_norm2);
		printf(" - L2 error norm:      % le\n", l2_err_norm2);

		if (h1_err_norm2 > EPS || l2_err_norm2 > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

#ifdef OUTPUT_DIR
		// output
		const char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			GmshOutputEngine output(ofile);
			output.out(&sln1, "Uh_1");
			output.out(&esln1, "U1");
			output.out(&sln2, "Uh_2");
			output.out(&esln2, "U2");

			fclose(ofile);
		}
		else {
			ERROR("Can not not open '%s' for writing.", of_name);
		}
#endif
	}

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

