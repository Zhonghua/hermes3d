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

//

EBCType bc_types(int marker) {
	// TODO:
	return BC_ESSENTIAL;
}

double bc_values(int marker, double x, double y, double z) {
	// TODO:
	return 0.0;
}

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	// TODO:
	return 0.0;
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	// TODO:
	return 0.0;
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

//	if (argc < 5) {
//		ERROR("Not enough parameters");
//		return ERR_NOT_ENOUGH_PARAMS;
//	}

//	sscanf(args[2], "%d", &m);
//	sscanf(args[3], "%d", &n);
//	sscanf(args[4], "%d", &o);

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
	space.set_bc_values(bc_values);

	// FIXME
//	int mx = maxn(4, m, n, o, 4);
//	order3_t order(mx, mx, mx);
//	printf("  - Setting uniform order to (%d, %d, %d)\n", mx, mx, mx);
//	space.set_uniform_order(order);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

#if defined WITH_SLEPC
	SlepcMatrix mat_a;
	SlepcMatrix mat_b;
	SlepcEigenSolver solver(mat_a, mat_b);
#endif

	// A
	Discretization da;
	da.set_num_equations(1);
	da.set_spaces(1, &space);
	da.set_pss(1, &pss);

	da.set_bilinear_form(0, 0, bilinear_form);
//	da.set_linear_form(0, linear_form);

	// assemble stiffness matrix
	da.create(&mat_a, NULL);

	Timer assemble_timer_a("Assembling stiffness matrix");
	assemble_timer_a.start();
	da.assemble(&mat_a, NULL);
	assemble_timer_a.stop();

	// B
	Discretization db;
	db.set_num_equations(1);
	db.set_spaces(1, &space);
	db.set_pss(1, &pss);

	db.set_bilinear_form(0, 0, bilinear_form);
//	db.set_linear_form(0, linear_form);

	// assemble stiffness matrix
	db.create(&mat_b, NULL);

	Timer assemble_timer_b("Assembling stiffness matrix");
	assemble_timer_b.start();
	db.assemble(&mat_b, NULL);
	assemble_timer_b.stop();


	// solve the stiffness matrix
	Timer solve_timer("Solving stiffness matrix");
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();

	// output the measured values
	printf("%s: %s (%lf secs)\n", assemble_timer_a.get_name(), assemble_timer_a.get_human_time(),
	    assemble_timer_a.get_seconds());
	printf("%s: %s (%lf secs)\n", assemble_timer_b.get_name(), assemble_timer_b.get_human_time(),
	    assemble_timer_b.get_seconds());
	printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(),
	    solve_timer.get_seconds());

	if (solved) {

#ifdef OUTPUT_DIR
#endif
	}
	else
		res = ERR_FAILURE;

#ifdef WITH_SLEPC
	mat_a.free();
	mat_b.free();
	SlepcFinalize();
#endif

	TRACE_END;

	return res;
}

