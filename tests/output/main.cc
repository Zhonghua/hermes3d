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

#include "config.h"
#include <math.h>
#include <hermes3d.h>

int m, n, o;

// error should be smaller than this epsilon
#define EPSILON								10e-10F

double fnc(double x, double y, double z) {
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

double dfnc(double x, double y, double z) {
	double ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	double ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	double ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);

	return -(ddxx + ddyy + ddzz);
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	// u(x, y, z) = pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4)
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fnc(x, y, z);
}

//

EBCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

double bc_values(int marker, double x, double y, double z) {
	return fnc(x, y, z);
}

scalar bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar linear_form(RealFunction *fv, RefMap *rv) {
	return int_F_v(dfnc, fv, rv);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	if (argc < 3) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	char *type = args[1];

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[2], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", args[2]);
		return ERR_FAILURE;
	}

#if defined GMSH
	GmshOutputEngine output(stdout);
#elif defined VTK
	VtkOutputEngine output(stdout);
#endif

	if (strcmp(type, "sln") == 0) {
		// Testing on Exact solution which always gives the same value (values from Solution may differ by epsilon)
		ExactSolution ex_sln(&mesh, exact_solution);

		output.out(&ex_sln, "U");
	}
	else if (strcmp(type, "ord") == 0) {
		H1ShapesetLobattoHex shapeset;
		PrecalcShapeset pss(&shapeset);

		H1Space space(&mesh, &shapeset);
		space.set_bc_types(bc_types);
		space.set_bc_values(bc_values);

		order3_t order(2, 3, 4);
		space.set_uniform_order(order);

		space.assign_dofs();

		output.out_orders(&space, "orders");
	}
	else if (strcmp(type, "bc") == 0) {
		output.out_bc(&mesh);
	}

	TRACE_END;

	return 0;
}

