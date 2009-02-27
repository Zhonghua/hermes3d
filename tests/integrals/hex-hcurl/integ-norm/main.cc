// This file is part of Hermes3D
//
// Copyright (c) 2009 Pavel Kus <pavel.kus@gmail.com>
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
 * main.cc
 *
 * Test of integrals and norms
 * Norm of a shape function is calculated using two different ways
 * via integrals and and norm
 * this test tests correctness of hcurl integrals (u.v, curl u . curl v),
 * hcurl and L2 norm in hcurl
 * and also transformations of solution and gradient in Solution
 * (because in each aproach, the transformation is done differently. When using
 * integrals, we transform curl using relation from book (Monk). When using norms,
 * partial derivatives are transformed in Solution and then curl is calculated in norm
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

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	//return int_grad_u_grad_v(fu, fv, ru, rv);
	//return 0;
	return rand();
}

double f(double x, double y, double z) {
	return 0;
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	//return int_F_v(f, fv, rv);
	return 0;
}

int main(int argc, char **args) {
	int res = ERR_SUCCESS;


#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

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

	int element_order = MAKE_HEX_ORDER(3, 4, 5);
	space.set_uniform_order(element_order);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

	UMFPackLinearSolver solver;

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
	Solution sln(&mesh);
	d.solve_system(1, &sln);

	scalar *sln_vector = new scalar[ndofs + 1];
	memset(sln_vector, 0, (ndofs + 1) * sizeof(scalar));

	int *fn_index = shapeset.get_bubble_indices(element_order);
	assert(shapeset.get_num_bubble_fns(element_order) == ndofs);

	Quad3D *quad = get_quadrature(MODE_HEXAHEDRON);
	RefMap rm(&mesh);
	rm.set_active_element(mesh.elements[0]);
	rm.set_quad(quad);

	for (int dof = 1; dof < ndofs + 1; dof++) {
		printf("processing dof %d ... ", dof);
		sln_vector[dof] = 1.;
		if (dof) sln_vector[dof - 1] = 0.;
		sln.set_solution_vector(sln_vector, false);

		pss.set_active_shape(fn_index[dof - 1]);

		scalar integ_l2_val = hcurl_int_u_v(&pss, &pss, &rm, &rm);
		scalar integ_hcurl_val = hcurl_int_curl_u_curl_v(&pss, &pss, &rm, &rm) + integ_l2_val;

		integ_l2_val = sqrt(integ_l2_val);
		integ_hcurl_val = sqrt(integ_hcurl_val);

		scalar norm_l2_val = hcurl_l2_norm(&sln);
		scalar norm_hcurl_val = hcurl_norm(&sln);

		if ((ABS(integ_l2_val - norm_l2_val) > EPS) || (ABS(integ_hcurl_val - norm_hcurl_val) > EPS)) {
			printf("ERROR l2 integ: " SPS ", norm: " SPS " hcurl integ: " SPS ", norm: " SPS "\n", SP(integ_l2_val), SP(norm_l2_val), SP(integ_hcurl_val), SP(norm_hcurl_val));
			res = ERR_FAILURE;
		}
		else
			printf("OK\n");
	}

	// destroy the solution vector
	delete[] sln_vector;

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

