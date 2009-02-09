/*
 * hex.cc
 *
 *
 *
 */
#include <math.h>

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

std::complex<double> imag(0, -1);

		///**********************************************************/
		///                      base 1
		///**********************************************************/

/*
double alpha = 1.0;

double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz, int comp) {
	switch(comp){
		case 0 : {
			dx = 0.;
			dy = -3./4. * y * (1. - z*z);
			dz = -3./4. * z * (1. - y*y);
			return 3./8. * (1. - y*y) * (1. - z*z);
		}
		case 1 : dx = dy = dz = 0; return 0;
		case 2 : dx = dy = dz = 0; return 0;
	}
}

scalar f(double x, double y, double z, int comp) {
	switch(comp){
		case 0 : return 3./4. * (2. - y*y - z*z) - alpha * (3./8. * (1. - y*y) * (1. - z*z));
		case 1 : return 0.;
		case 2 : return 0.;
	}
}

double bc_values(int marker, double x, double y, double z, int comp) {
	switch (marker) {
		case 0: return 0.;
		case 1: {
			switch(comp){
				case 0 : return 0.;
				case 1 : return 3./4. * y * (1 - z*z);
				case 2 : return 3./4. * z * (1 - y*y);
			}
		}
		case 2: return 0.;
		case 3: return 0.;
		case 4: return 0.;
		case 5: return 0.;
		default: EXIT(ERR_FACE_INDEX_OUT_OF_RANGE);
	}
}

EBCType bc_types(int marker) {
//	if(marker == 1)
//		return BC_NATURAL;
//	else
		return BC_ESSENTIAL;
}
*/

		///**********************************************************/
		///              base 2 multiplied by i
		///**********************************************************/

std::complex<double> alpha(-1, 1);

scalar exact_solution(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz, int comp) {
	switch(comp){
		case 0 : {
			dx = imag * 3./8. * (1. - y*y) * (1. - z*z);
			dy = imag * (-3./4.) * x * y * (1. - z*z);
			dz = imag * (-3./4.) * x * (1. - y*y) * z;
			return  imag * 3./8. * x * (1. - y*y) * (1. - z*z);
		}
		case 1 :return 0;
		case 2 :return 0;
	}
}

scalar f(double x, double y, double z, int comp) {
	switch(comp){
		case 0 : return imag * (3./4. * x * (2. - y*y - z*z) - alpha * (3./8. * x * (1. - y*y) * (1. - z*z)));
		case 1 : return imag * (-3./4. * y * (1. - z*z));
		case 2 : return imag * (-3./4. * z * (1. - y*y));
	}
}

double bc_values(int marker, double x, double y, double z, int comp) {
	return 0.;
}


EBCType bc_types(int marker) {
		return BC_ESSENTIAL;
}


		///**********************************************************/
		///            general polynomial function
      ///            satisfying perfect condurctor bc
		///**********************************************************/

/*
double alpha = 1.0;

double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz, int comp) {
	switch(comp){
		case 0 : {
			dx = 0.;
			dy = -2*y*(1-z*z);
			dz = -2*z*(1-y*y);
			return (1-y*y) * (1-z*z);
		}
		case 1 : {
			dx = -2*x*(1-z*z);
			dy = 0.;
			dz = -2*z*(1-x*x);
			return (1-x*x) * (1-z*z);
		}
		case 2 : {
			dx = -2*x*(1-y*y);
			dy = -2*y*(1-x*x);
			dz = 0.;
			return (1-x*x) * (1-y*y);
		}
	}
}

scalar f(double x, double y, double z, int comp) {
	scalar curlpart;
	double dx, dy, dz;
	switch(comp){
		case 0 : curlpart = 4 - 2*y*y - 2*z*z; break;
		case 1 : curlpart = 4 - 2*x*x - 2*z*z; break;
		case 2 : curlpart = 4 - 2*x*x - 2*y*y; break;
	}

	return curlpart - alpha * exact_solution(x, y, z, dx, dy, dz, comp);
}

double bc_values(int marker, double x, double y, double z, int comp) {
	return 0;
}

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}
*/
		///**********************************************************/
		///            general polynomial function
      ///            satisfying perfect condurctor bc
		///**********************************************************/

/*
std::complex<double> alpha(-1, 1);

scalar exact_solution(double x, double y, double z, scalar &dx, scalar &dy, scalar &dz, int comp) {
	switch(comp){
		case 0 : {
			dx = 3.*x*x;
			dy = 0;
			dz = imag * 3.*z*z;
			return x*x*x + imag * z*z*z;
		}
		case 1 : {
			dx = y*z + imag * 2.*x*y*y;
			dy = x*z + imag * 2.*x*x*y;
			dz = x*y;
			return x*y*z + imag * x*x*y*y;
		}
		case 2 : {
			dx = -imag * y*y*y*z;
			dy = -imag * 3.*x*y*y*z;
			dz = 1. - imag * x*y*y*y;
			return z - imag * x*y*y*y*z;
		}
	}
}

scalar f(double x, double y, double z, int comp) {
	scalar curlpart;
	scalar dx, dy, dz;
	switch(comp){
		case 0 : curlpart = z + imag * (4.*x*y - 6.*z + y*y*y); break;
		case 1 : curlpart = -imag * y*y * (3.*x + 2.); break;
		case 2 : curlpart = x + imag * 6.*x*y*z; break;
	}

	return curlpart - alpha * exact_solution(x, y, z, dx, dy, dz, comp);
}

double bc_values(int marker, double x, double y, double z, int comp) {
	return 0;
}

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}*/

	///**********************************************************/
	///              definition of the forms
	///**********************************************************/

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
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

