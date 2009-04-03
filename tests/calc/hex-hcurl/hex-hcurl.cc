// This file is part of Hermes3D
//
// Copyright (c) 2008 - 2009 Pavel Kus <pavel.kus@gmail.com>
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

// general polynomial function satisfying perfect conductor bc
// exact solution has zero tangential component on the boundary on the domain (-1, 1)^3
// DO NOT TEST on other domains than (-1, 1)^3

const double alpha = 1.0;

scalar3 &exact_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz) {
	dx[0] = (1 - y*y)*(1 - z*z)*(z - 6*x*x);
	dx[1] = 2*(1 - x*x)*(1 - z*z) - 2*x*(1 - z*z)*(y*y*y + 2*x);
	dx[2] = -2*x*(1 - y*y)*(z*z - 3*x*y*z) - 3*y*z*(1 - x*x)*(1 - y*y);

	dy[0] = -2*y*(1 - z*z)*(1 - 2*x*x*x + x*z);
	dy[1] = 3*y*y*(1 - x*x)*(1 - z*z);
	dy[2] = -2*y*(1 - x*x)*(z*z - 3*x*y*z) - 3*x*z*(1 - x*x)*(1 - y*y);

	dz[0] = x*(1 - y*y)*(1 - z*z) - 2*z*(1 - y*y)*(1 - 2*x*x*x + x*z);
	dz[1] = -2*z*(1 - x*x)*(y*y*y + 2*x);
	dz[2] = (1 - x*x)*(1 - y*y)*(2*z - 3*x*y);

	static scalar3 val;
	val[0] = (1-y*y) * (1-z*z) * (x*z - 2*x*x*x + 1);
	val[1] = (1-x*x) * (1-z*z) * (y*y*y + 2*x);
	val[2] = (1-x*x) * (1-y*y) * (z*z - 3*x*y*z);

	return val;
}

scalar3 &f(double x, double y, double z) {
	scalar3 curlpart = {
		2*(1 - y*y)*(1 - 2*x*x*x + x*z) + 2*(1 - z*z)*(1 - 2*x*x*x + x*z) - 6*x*y*y*(1 - z*z) - 3*y*(1 - x*x)*(1 - y*y) - 2*x*(1 - y*y)*(2*z - 3*x*y) + 4*x*z*(1 - y*y),
		2*(1 - x*x)*(y*y*y + 2*x) + 2*(1 - z*z)*(y*y*y + 2*x) + 8*x*(1 - z*z) - 3*x*(1 - x*x)*(1 - y*y) - 2*y*(1 - x*x)*(2*z - 3*x*y) - 2*y*(1 - z*z)*(z - 6*x*x),
		(1 - y*y)*(1 - z*z) + 2*(1 - x*x)*(z*z - 3*x*y*z) + 2*(1 - y*y)*(z*z - 3*x*y*z) - 6*z*y*y*(1 - x*x) - 2*z*(1 - y*y)*(z - 6*x*x) - 12*x*y*z*(1 - x*x) - 12*x*y*z*(1 - y*y)
	};

	scalar3 dx, dy, dz;
	scalar3 &ev = exact_solution(x, y, z, dx, dy, dz);

	static scalar3 val;
	val[0] = curlpart[0] - alpha * ev[0];
	val[1] = curlpart[1] - alpha * ev[1];
	val[2] = curlpart[2] - alpha * ev[2];

	return val;
}

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

// definition of the forms

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return hcurl_int_curl_u_curl_v(fu, fv, ru, rv) - alpha * hcurl_int_u_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	return hcurl_int_F_v(f, fv, rv);
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

	if (argc < 3) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	HcurlShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[1]);
		return ERR_FAILURE;
	}

	printf("* Setting the space up\n");
	HcurlSpace space(&mesh, &shapeset);
	space.set_bc_types(bc_types);

	int o;
	sscanf(args[2], "%d", &o);
	order3_t order(o, o, o);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o, o, o);
	space.set_uniform_order(order);

	int ndofs = space.assign_dofs();
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
	d.set_num_equations(1);
	d.set_spaces(1, &space);
	d.set_pss(1, &pss);

	d.set_bilinear_form(0, 0, bilinear_form);
	d.set_linear_form(0, linear_form);

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

	if (solved) {
		Solution sln(&mesh);
		sln.set_space_and_pss(&space, &pss);
		sln.set_solution_vector(solver.get_solution(), false);

		printf("* Solution:\n");
		scalar *s = sln.get_solution_vector();
		for (int i = 1; i <= ndofs; i++) {
			printf(" x[% 3d] = " SCALAR_FMT "\n", i, SCALAR(s[i]));
		}

		// output the measured values
		printf("%s: %s (%lf secs)\n", assemble_timer.get_name(), assemble_timer.get_human_time(), assemble_timer.get_seconds());
		printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(), solve_timer.get_seconds());

		ExactSolution ex_sln(&mesh, exact_solution);

		double hcurl_sln_norm = hcurl_norm(&sln);
		double hcurl_err_norm = hcurl_error(&sln, &ex_sln);
		printf(" - Hcurl solution norm:   % le\n", hcurl_sln_norm);
		printf(" - Hcurl error norm:      % le\n", hcurl_err_norm);

		double l2_sln_norm = l2_norm_hcurl(&sln);
		double l2_err_norm = l2_error_hcurl(&sln, &ex_sln);
		printf(" - L2 solution norm:   % le\n", l2_sln_norm);
		printf(" - L2 error norm:      % le\n", l2_err_norm);

		if (hcurl_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}


#ifdef OUTPUT_DIR
		// output
		printf("starting output\n");
		const char *of_name = OUTPUT_DIR "/solution.vtk";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			ExactSolution ex_sln(&mesh, exact_solution);
			DiffFilter eh(&sln, &ex_sln);
			DiffFilter eh_dx(&sln, &ex_sln, FN_DX, FN_DX);
//			DiffFilter eh_dy(mesh, &sln, &ex_sln, FN_DY, FN_DY);
//			DiffFilter eh_dz(mesh, &sln, &ex_sln, FN_DZ, FN_DZ);

//			GmshOutputEngine output(ofile);
			VtkOutputEngine output(ofile);
			output.out(&sln, "Uh", FN_VAL);
//			output.out(&sln, "Uh_0", FN_VAL_0);
//			output.out(&sln, "Uh_1", FN_VAL_1);
//			output.out(&sln, "Uh_2", FN_VAL_2);
//			output.vector_out(&sln, "Uh_dx", FN_DX);
//			output.vector_out(&sln, "Uh_dy", FN_DY);
//			output.vector_out(&sln, "Uh_dz", FN_DZ);
//			output.scalar_out(&sln, "Uh_dx_0", FN_DX_0);
//			output.scalar_out(&sln, "Uh_dx_1", FN_DX_1);
//			output.scalar_out(&sln, "Uh_dx_2", FN_DX_2);
//			output.out(&sln, "Uh dy", FN_DY_0);
//			output.out(&sln, "Uh dz", FN_DZ_0);

//			output.vector_out(&sln, "Uh_dx", FN_DX);
//			output.vector_out(&sln, "Uh_dy", FN_DY);
//			output.vector_out(&sln, "Uh_dz", FN_DZ);
//			output.scalar_out(&sln, "Uh_dx_0", FN_DX_0);
//			output.scalar_out(&sln, "Uh_dx_1", FN_DX_1);
//			output.scalar_out(&sln, "Uh_dx_2", FN_DX_2);
//			output.out(&sln, "Uh_dy", FN_DY_0);
//			output.out(&sln, "Uh_dz", FN_DZ_0);
//			output.out(&eh, "Eh", FN_VAL);
//			output.out(&eh_dx, "Eh_dx", FN_VAL);
//			output.out(&eh_dy, "Eh_dy");
//			output.out(&eh_dz, "Eh_dz");
//			output.out(&ex_sln, "U", FN_VAL);
//			output.out(&ex_sln, "U_dx", FN_DX);
//			output.out(&ex_sln, "U_dy", FN_DY);
//			output.out(&ex_sln, "U_dz", FN_DZ);
//			output.scalar_out(&ex_sln, "U_0", FN_VAL_0);
//			output.scalar_out(&ex_sln, "U_1", FN_VAL_1);
//			output.scalar_out(&ex_sln, "U_2", FN_VAL_2);
//			output.out(&ex_sln, "U_dy", FN_DY_0);
//			output.out(&ex_sln, "U_dz", FN_DZ_0);
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

