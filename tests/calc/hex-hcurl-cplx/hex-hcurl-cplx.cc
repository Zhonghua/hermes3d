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

#include <math.h>

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

std::complex<double> imag(0, -1);
std::complex<double> alpha(-1, 1);

scalar3 &exact_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz) {
	dx[0] = imag * 3./8. * (1. - y*y) * (1. - z*z);
	dx[1] = 0.0;
	dx[2] = 0.0;

	dy[0] = imag * (-3./4.) * x * y * (1. - z*z);
	dy[1] = 0.0;
	dy[2] = 0.0;

	dz[0] = imag * (-3./4.) * x * (1. - y*y) * z;
	dz[1] = 0.0;
	dz[2] = 0.0;

	static scalar3 val;
	val[0] = imag * 3./8. * x * (1. - y*y) * (1. - z*z);
	val[1] = 0.0;
	val[2] = 0.0;
	return val;
}

scalar3 &f(double x, double y, double z) {
	static scalar3 val;
	val[0] = imag * (3./4. * x * (2. - y*y - z*z) - alpha * (3./8. * x * (1. - y*y) * (1. - z*z)));
	val[1] = imag * (-3./4. * y * (1. - z*z));
	val[2] = imag * (-3./4. * z * (1. - y*y));
	return val;
}

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

/// definition of the forms

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return hcurl_int_curl_u_curl_v(fu, fv, ru, rv) - alpha * hcurl_int_u_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	return hcurl_int_F_v(f, fv, rv);
}

scalar bilinear_form_surf(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp) {
	return -hcurl_surf_int_u_v(fu, fv, ru, rv, fp);
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

	HcurlShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", argv[1]);
		return ERR_FAILURE;
	}

	printf("* Setting the space up\n");
	HcurlSpace space(&mesh, &shapeset);
	space.set_bc_types(bc_types);

	int order;
	sscanf(argv[2], "%d", &order);
	int dir_x = order, dir_y = order, dir_z = order;
	order3_t o(dir_x, dir_y, dir_z);
	printf("  - Setting uniform order to (%d, %d, %d)\n", dir_x, dir_y, dir_z);
	space.set_uniform_order(o);

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

	d.set_bilinear_form(0, 0, bilinear_form, NULL, bilinear_form_surf);
	d.set_linear_form(0, linear_form, linear_form_surf);

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

#ifdef OUTPUT_DIR
	solver.dump_matrix(OUTPUT_DIR "/matrix");
#endif

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

		// norm
		ExactSolution ex_sln(&mesh, exact_solution);
		double hcurl_sln_norm = hcurl_norm(&sln);
		double hcurl_err_norm = hcurl_error(&sln, &ex_sln);
		printf(" - Hcurl solution norm: % le\n", hcurl_sln_norm);
		printf(" - Hcurl error norm:    % le\n", hcurl_err_norm);

		double l2_sln_norm = l2_norm_hcurl(&sln);
		double l2_err_norm = l2_error_hcurl(&sln, &ex_sln);
		printf(" - L2 solution norm:    % le\n", l2_sln_norm);
		printf(" - L2 error norm:       % le\n", l2_err_norm);

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
			ExactSolution ex_sln(&mesh, exact_solution_0, exact_solution_1, exact_solution_2);

			RealPartFilter real_sln(&mesh, &sln, FN_VAL);
			ImagPartFilter imag_sln(&mesh, &sln, FN_VAL);

			DiffFilter eh(&mesh, &sln, &ex_sln);
			DiffFilter eh_dx(&mesh, &sln, &ex_sln, FN_DX, FN_DX);
//			DiffFilter eh_dy(&mesh, &sln, &ex_sln, FN_DY, FN_DY);
//			DiffFilter eh_dz(&mesh, &sln, &ex_sln, FN_DZ, FN_DZ);

//			GmshOutputEngine output(ofile);
			VtkOutputEngine output(ofile);

			output.out(&real_sln, "real_Uh", FN_VAL);
			output.out(&imag_sln, "imag_Uh", FN_VAL);

			output.out(&real_sln, "real_Uh_0", FN_VAL_0);
			output.out(&real_sln, "real_Uh_1", FN_VAL_1);
			output.out(&real_sln, "real_Uh_2", FN_VAL_2);

			output.out(&imag_sln, "imag_Uh_0", FN_VAL_0);
			output.out(&imag_sln, "imag_Uh_1", FN_VAL_1);
			output.out(&imag_sln, "imag_Uh_2", FN_VAL_2);

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

