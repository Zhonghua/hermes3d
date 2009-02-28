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

/*
 * hex-hcurl-imped.cc
 * Time harmonic Maxwell equations with impedance BC
 *
 * calculates solution of the problem:
 *  curl(curl(u)) = f on Omega
 *  u x n = 0 on Gamma_P
 *  curl(u) x n - i (n x u) x n = g on Gamma_I
 *  where n is outer normal
 */

#include "config.h"
#ifdef WITH_PETSC
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

std::complex<double> imag(0, 1);

		///**********************************************************/
		///                    E = (1, 0, 0)
		///**********************************************************/
/*
scalar exact_solution(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz, int comp) {
	switch(comp){
		case 0 : dx = dy = dz = 0.; return 1.;
		case 1 : dx = dy = dz = 0.; return 0.;
		case 2 : dx = dy = dz = 0.; return 0.;
}
}

scalar f(double x, double y, double z, int comp) {
	switch(comp){
		case 0 : return -1.;
		case 1 : return 0.;
		case 2 : return 0.;
}
}

scalar bc_values(int marker, double x, double y, double z, int comp) {
	scalar bc[3] = {0., 0., 0.};
	switch (marker) {
		case 0: break;
		case 1: break;
		case 2: bc[0] = -imag; break;
		case 3: bc[0] = -imag; break;
		case 4: bc[0] = -imag; break;
		case 5: bc[0] = -imag; break;
		default:assert(0); //EXIT(ERR_FACE_INDEX_OUT_OF_RANGE);
	}
	return bc[comp];
}

EBCType bc_types(int marker) {
	if((marker == 2) || (marker == 3) || (marker == 4) || (marker == 5))
		return BC_NATURAL;
	else
		return BC_ESSENTIAL;
}
*/
		///**********************************************************/
		///                    E = (0, x, 0)
		///**********************************************************/
/*
scalar exact_solution(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz, int comp) {
	switch(comp){
		case 0 : dx = dy = dz = 0; return 0;
		case 1 : dx = 1.; dy = dz = 0; return x;
		case 2 : dx = dy = dz = 0; return 0;
	}
}

scalar f(double x, double y, double z, int comp) {
	switch(comp){
		case 0 : return 0.;
		case 1 : return -x;
		case 2 : return 0.;
	}
}

scalar bc_values(int marker, double x, double y, double z, int comp) {
	scalar bc[3] = {0., 0., 0.};
	switch (marker) {
		case 0: bc[1] = -1. - imag * x; break;
		case 1: bc[1] =  1. - imag * x; break;
		case 2: break;
		case 3: break;
		case 4: bc[1] = - imag * x; break;
		case 5: bc[1] = - imag * x; break;
		default: EXIT(ERR_FACE_INDEX_OUT_OF_RANGE);
	}
	return bc[comp];
}

EBCType bc_types(int marker) {
	if((marker == 0) || (marker == 1) || (marker == 4) || (marker == 5))
		return BC_NATURAL;
	else
		return BC_ESSENTIAL;
}
*/
		///**********************************************************/
		///                    E = (0, x*x, 0)
		///**********************************************************/
/*
scalar exact_solution(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz, int comp) {
	switch(comp){
		case 0 : dx = dy = dz = 0; return 0;
		case 1 : dx = 2*x; dy = dz = 0; return x*x;
		case 2 : dx = dy = dz = 0; return 0;
	}
}

scalar f(double x, double y, double z, int comp) {
	switch(comp){
		case 0 : return 0.;
		case 1 : return -2. - x*x;
		case 2 : return 0.;
	}
}

scalar bc_values(int marker, double x, double y, double z, int comp) {
	scalar bc[3] = {0., 0., 0.};
	switch (marker) {
		case 0: bc[1] = (-2.*x) - imag * x*x; break;
		case 1: bc[1] = (2.*x)  - imag * x*x; break;
		case 2: break;
		case 3: break;
		case 4: bc[1] = - imag * x*x; break;
		case 5: bc[1] = - imag * x*x; break;
		default: EXIT(ERR_FACE_INDEX_OUT_OF_RANGE);
	}
	return bc[comp];
}

EBCType bc_types(int marker) {
	if((marker == 0) || (marker == 1) || (marker == 4) || (marker == 5))
		return BC_NATURAL;
	else
		return BC_ESSENTIAL;
}
*/
		///**********************************************************/
		///            general polynomial function
		///**********************************************************/
scalar exact_solution(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz, int comp) {
	switch(comp){
		case 0 : {
			dx = 3.*x*x*y*y - 3.*y*y*y*z;
			dy = 2.*x*x*x*y - 9.*x*y*y*z;
			dz = -3.*x*y*y*y;
			return x*x*x*y*y - 3.*x*y*y*y*z;
		}
		case 1 : {
			dx = 3.*x*x*y*y*y*z*z*z + 4.*x*y*z;
			dy = 3.*x*x*x*y*y*z*z*z + 2.*x*x*z;
			dz = 3.*x*x*x*y*y*y*z*z + 2.*x*x*y;
			return x*x*x*y*y*y*z*z*z + 2.*x*x*y*z;
		}
		case 2 : {
			dx = -12.*x*x;
			dy = z*z*z;
			dz = 3.*y*z*z;
			return y*z*z*z - 4.*x*x*x;
		}
	}
}

scalar f(double x, double y, double z, int comp) {
	scalar curlpart;
	scalar dx, dy, dz;
	switch(comp){
		case 0 : curlpart = 4*x*z + 18*x*y*z - 2*x*x*x + 9*x*x*y*y*z*z*z; break;
		case 1 : curlpart = -4*y*z + 3*z*z - 9*z*y*y + 6*y*x*x - 6*x*y*y*y*z*z*z - 6*z*x*x*x*y*y*y; break;
		case 2 : curlpart = 24*x + 2*x*x + 9*x*x*x*y*y*z*z - 3*y*y*y; break;
	}

	return curlpart - exact_solution(x, y, z, dx, dy, dz, comp);
}

//TODO this could be writen in much simplier way. Just use curl of exact solution
//and cross product defined in Scalar3D..
scalar bc_values(int marker, double x, double y, double z, int comp) {
	scalar bc[3] = {0., 0., 0.};
	switch (marker) {
		case 0:
			bc[1] = -4*x*y*z - 9*x*z*y*y + 2*y*x*x*x - 3*x*x*y*y*y*z*z*z;
			bc[2] = 12*x*x - 3*x*y*y*y;
			break;
		case 1:
			bc[1] = 4*x*y*z + 9*x*z*y*y - 2*y*x*x*x + 3*x*x*y*y*y*z*z*z;
			bc[2] = -12*x*x + 3*x*y*y*y;
			break;
		case 2:
			bc[0] = 4*x*y*z + 9*x*z*y*y - 2*y*x*x*x + 3*x*x*y*y*y*z*z*z;
			bc[2] = 2*y*x*x + 3*x*x*x*y*y*y*z*z - z*z*z;
			break;
		case 3:
			bc[0] = -4*x*y*z - 9*x*z*y*y + 2*y*x*x*x - 3*x*x*y*y*y*z*z*z;
			bc[2] = -2*y*x*x - 3*x*x*x*y*y*y*z*z + z*z*z;
		break;
		case 4:
			bc[0] = -12*x*x + 3*x*y*y*y;
			bc[1] = -2*y*x*x - 3*x*x*x*y*y*y*z*z + z*z*z;
			break;
		case 5:
			bc[0] = 12*x*x - 3*x*y*y*y;
			bc[1] = 2*y*x*x + 3*x*x*x*y*y*y*z*z - z*z*z;
			break;
		default: EXIT(ERR_FACE_INDEX_OUT_OF_RANGE);
	}
	switch (marker){
		case 0:
		case 1:
			bc[1] -= imag * (2*y*z*x*x + x*x*x*y*y*y*z*z*z);
			bc[2] -= imag * (-4*x*x*x + y*z*z*z);
			break;
		case 2:
		case 3:
			bc[0] -= imag * (x*x*x*y*y - 3*x*z*y*y*y);
			bc[2] -= imag * (-4*x*x*x + y*z*z*z);
			break;
		case 4:
		case 5:
			bc[0] -= imag * (x*x*x*y*y - 3*x*z*y*y*y);
			bc[1] -= imag * (2*y*z*x*x + x*x*x*y*y*y*z*z*z);
			break;
		default: EXIT(ERR_FACE_INDEX_OUT_OF_RANGE);
	}
	return bc[comp];
}

EBCType bc_types(int marker) {
	return BC_NATURAL;
}



	///**********************************************************/
	///              definition of the forms
	///**********************************************************/

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return hcurl_int_curl_u_curl_v(fu, fv, ru, rv) - hcurl_int_u_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	return hcurl_int_F_v(f, fv, rv);
}

scalar bilinear_form_surf(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp) {
	return -imag * hcurl_surf_int_u_v(fu, fv, ru, rv, fp);
}

scalar linear_form_surf(RealFunction *fv, RefMap *rv, FacePos *fp) {
	return hcurl_surf_int_G_v(fv, rv, fp);
}

scalar exact_solution_0(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz) {
	return exact_solution(x, y, z, dx, dy, dz, 0);
}

scalar exact_solution_1(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz) {
	return exact_solution(x, y, z, dx, dy, dz, 1);
}

scalar exact_solution_2(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz) {
	return exact_solution(x, y, z, dx, dy, dz, 2);
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

	HCurlShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[1]);
		return ERR_FAILURE;
	}

	printf("* Setting the space up\n");
	HCurlSpace space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_bc_values(bc_values);

	int order;
	sscanf(args[2], "%d", &order);
	int dir_x = order, dir_y = order, dir_z = order;
	int o = MAKE_HEX_ORDER(dir_x, dir_y, dir_z);
	printf("  - Setting uniform order to (%d, %d, %d)\n", dir_x, dir_y, dir_z);
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

	d.set_bilinear_form(0, 0, bilinear_form, NULL, bilinear_form_surf);
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

	solver.dump_matrix("output/matrix");
	solver.dump_rhs("output/rhs");

	if (solved) {
		printf("* Solution:\n");
		scalar *s = sln.get_solution_vector();
		for (int i = 1; i <= ndofs; i++) {
			printf(" x[% 3d] = " SPS "\n", i, SP(s[i]));
		}

		// output the measured values
		printf("%s: %s (%lf secs)\n", assemble_timer.get_name(), assemble_timer.get_human_time(), assemble_timer.get_seconds());
		printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(), solve_timer.get_seconds());

		// norm
		double hcurl_sln_norm = hcurl_norm(&sln);
		double hcurl_err_norm = hcurl_error_norm_exact(&sln, exact_solution);
		printf(" - HCurl solution norm:   % le\n", hcurl_sln_norm);
		printf(" - HCurl error norm:      % le\n", hcurl_err_norm);

		double l2_sln_norm = hcurl_l2_norm(&sln);
		double l2_err_norm = hcurl_l2_error_norm_exact(&sln, exact_solution);
		printf(" - L2 solution norm:   % le\n", l2_sln_norm);
		printf(" - L2 error norm:      % le\n", l2_err_norm);

		if (hcurl_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}


#ifdef OUTPUT_DIR
		// output
		printf("starting output\n");
		const char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			ExactSolution ex_sln(&mesh, exact_solution_0, exact_solution_1, exact_solution_2);

			RealPartFilter real_sln(&mesh, &sln, FN_VAL);
			ImagPartFilter imag_sln(&mesh, &sln, FN_VAL);

			GmshOutputEngine output(ofile);

			output.out(&real_sln, "real_Uh", FN_VAL);
			output.out(&imag_sln, "imag_Uh", FN_VAL);

			output.out(&real_sln, "real_Uh_0", FN_VAL_0);
			output.out(&real_sln, "real_Uh_1", FN_VAL_1);
			output.out(&real_sln, "real_Uh_2", FN_VAL_2);

			output.out(&imag_sln, "imag_Uh_0", FN_VAL_0);
			output.out(&imag_sln, "imag_Uh_1", FN_VAL_1);
			output.out(&imag_sln, "imag_Uh_2", FN_VAL_2);

			DiffFilter eh(&mesh, &sln, &ex_sln);
//			DiffFilter eh_dx(mesh, &sln, &ex_sln, FN_DX, FN_DX);
//			DiffFilter eh_dy(mesh, &sln, &ex_sln, FN_DY, FN_DY);
//			DiffFilter eh_dz(mesh, &sln, &ex_sln, FN_DZ, FN_DZ);

			RealPartFilter real_eh(&mesh, &eh, FN_VAL);
			ImagPartFilter imag_eh(&mesh, &eh, FN_VAL);

			output.out(&real_eh, "real_Eh", FN_VAL);
			output.out(&imag_eh, "imag_Eh", FN_VAL);

			fclose(ofile);
		}
		else {
			ERROR("Can not not open '%s' for writing.", of_name);
		}
#endif
	}

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

