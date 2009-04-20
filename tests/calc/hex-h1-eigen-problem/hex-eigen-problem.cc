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

#include "config.h"
#include <math.h>
#ifdef WITH_SLEPC
#include <slepc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>
#include <common/utils.h>

#define BEGIN_BLOCK							{
#define END_BLOCK							}

//

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

double bc_values(int marker, double x, double y, double z) {
	return 0.0;
}

scalar stiff_bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar mass_bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_u_v(fu, fv, ru, rv);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;

#ifdef WITH_SLEPC
	SlepcInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	if (argc < 3) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	int om;
	sscanf(args[2], "%d", &om);

	H1ShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[1]);
		return ERR_FAILURE;
	}

//	mesh.refine_all_elements(REFT_HEX_XYZ);
//	mesh.refine_all_elements(REFT_HEX_XYZ);
//	mesh.refine_all_elements(REFT_HEX_XYZ);

	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_bc_values(bc_values);

	order3_t order(om, om, om);
	space.set_uniform_order(order);

	int ndofs = space.assign_dofs();

	BEGIN_BLOCK
#if defined WITH_SLEPC
	PetscMatrix mat_a;
	PetscMatrix mat_b;
#endif

	// A (stiffness matrix)
	Discretization da;
	da.set_num_equations(1);
	da.set_spaces(1, &space);
	da.set_pss(1, &pss);

	da.set_bilinear_form(0, 0, stiff_bilinear_form);

	da.create(&mat_a, NULL);
	da.assemble(&mat_a, NULL);
	mat_a.finish();

	// B (mass matrix)
	Discretization db;
	db.set_num_equations(1);
	db.set_spaces(1, &space);
	db.set_pss(1, &pss);

	db.set_bilinear_form(0, 0, mass_bilinear_form);

	db.create(&mat_b, NULL);
	db.assemble(&mat_b, NULL);
	mat_b.finish();

	// solve the problem
#if defined WITH_SLEPC
	SlepcEigenSolver solver(mat_a, mat_b);
#endif

	Timer solve_timer("Solving the problem");
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();

	if (solved) {
		int nc = solver.get_converged();
		for (int i = 0; i < nc; i++) {
			printf("%d-th pair\n", i);

			double kr, ki;							// real and imaginary part of the eigenvalue
			double *xr = new double [ndofs + 1];	// real and imaginary part of the eigenvector
			double *xi = new double [ndofs + 1];
			solver.get_eigen_pair(i, &kr, &ki, xr, xi);

			printf("  Eigenvalue  : % lf + % lfi\n", kr, ki);
			printf("  Eigenvector : ");
			for (int j = 1; j <= ndofs; j++) {
				if (j > 1) printf(", ");
				printf("% lf", xr[j]);
			}
			printf("\n");

			double err = solver.compute_relative_error(i);
			printf("  Rel. error  : % e\n", err);

			Solution sln(&mesh);
			sln.set_space_and_pss(&space, &pss);
			sln.set_solution_vector(xr, false);

#ifdef OUTPUT_DIR
			// output
			char of_name[1024];
			sprintf(of_name, "%s/eigen-fn-%d.pos", OUTPUT_DIR, i);
			FILE *ofile = fopen(of_name, "w");
			if (ofile != NULL) {
				GmshOutputEngine output(ofile);
				output.out(&sln, "Uh", FN_VAL_0);

				fclose(ofile);
			}
			else {
				ERROR("Can't open '%s' for writing.", of_name);
			}
#endif

			delete [] xr;
			delete [] xi;
		}
	}
	else
		res = ERR_FAILURE;

	END_BLOCK
#ifdef WITH_SLEPC
	SlepcFinalize();
#endif

	TRACE_END;

	return res;
}

