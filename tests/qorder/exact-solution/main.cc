// This file is part of Hermes3D
//
// Copyright (c) 2008 - 2009 David Andrs <dandrs@unr.edu>
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

//
// main.cc
//
// Testing if the values calculated by ExactSolution class are correct.
// The purpose is to check if the Qorder stuff works okey.
//
//
// TODO:
// - more functions
// - test on vector-valued shapesets
//

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

// error should be smaller than this epsilon
#define EPS								10e-14F

double (*exact_solution)(double x, double y, double z, double &dx, double &dy, double &dz);


// FN #1

int m = 2, n = 2, o = 2;

double fn1_fnc(double x, double y, double z) {
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

double fn1_exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	// u(x, y, z) = pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4)
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fn1_fnc(x, y, z);
}

////////////////////////

bool test_vertex_values(MeshFunction &fu, RefMap &ru, Quad3D *&quad) {
	fu.set_quad_order(VTX_QORDER());

	double *val = fu.get_fn_values();
	int np = quad->get_vertex_num_points();

	double *x = ru.get_vertex_phys_x();
	double *y = ru.get_vertex_phys_y();
	double *z = ru.get_vertex_phys_z();

	for (int k = 0; k < np; k++) {
		double dx, dy, dz;
		double fn = exact_solution(x[k], y[k], z[k], dx, dy, dz);

		if (fabs(fn - val[k]) != 0.0) return false;
	}

	return true;
}

bool test_edge_values(int edge, order1_t order, MeshFunction &fu, RefMap &ru, Quad3D *&quad) {
	fu.set_quad_order(EDGE_QORDER(edge, order));

	double *val = fu.get_fn_values();
	int np = quad->get_edge_num_points(order);

	double *x = ru.get_edge_phys_x(edge, order);
	double *y = ru.get_edge_phys_y(edge, order);
	double *z = ru.get_edge_phys_z(edge, order);

	for (int k = 0; k < np; k++) {
		double dx, dy, dz;
		double fn = exact_solution(x[k], y[k], z[k], dx, dy, dz);

		if (fabs(fn - val[k]) != 0.0) return false;
	}

	return true;
}

bool test_face_values(int face, order2_t order, MeshFunction &fu, RefMap &ru, Quad3D *&quad) {
	fu.set_quad_order(FACE_QORDER(face, order));

	double *val = fu.get_fn_values();
	int np = quad->get_face_num_points(face, order);

	double *x = ru.get_face_phys_x(face, order);
	double *y = ru.get_face_phys_y(face, order);
	double *z = ru.get_face_phys_z(face, order);

	for (int k = 0; k < np; k++) {
		double dx, dy, dz;
		double fn = exact_solution(x[k], y[k], z[k], dx, dy, dz);

		if (fabs(fn - val[k]) != 0.0) return false;
	}

	return true;
}

bool test_elem_values(order3_t order, MeshFunction &fu, RefMap &ru, Quad3D *&quad) {
	fu.set_quad_order(ELEM_QORDER(order));

	double *val = fu.get_fn_values();
	int np = quad->get_num_points(order);

	double *x = ru.get_phys_x(order);
	double *y = ru.get_phys_y(order);
	double *z = ru.get_phys_z(order);

	for (int k = 0; k < np; k++) {
		double dx, dy, dz;
		double fn = exact_solution(x[k], y[k], z[k], dx, dy, dz);

		if (fabs(fn - val[k]) != 0.0) return false;
	}

	return true;
}

/////////

int main(int argc, char *argv[]) {
	int res = ERROR_SUCCESS;

	H1ShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", argv[1]);
		return ERR_FAILURE;
	}

	// fn1
	exact_solution = fn1_exact_solution;

	bool passed = true;
	for (int o = 1; o < MAX_QUAD_ORDER; o++) {
		order3_t order(o, o, o);
		printf("Order #%s\n", order.str());
		// test exact solution
		ExactSolution ex_sln(&mesh, fn1_exact_solution);
		FOR_ALL_ACTIVE_ELEMENTS(idx, &mesh) {
			Element *e = mesh.elements[idx];

			// use independent ref. map
			RefMap rm(&mesh);
			rm.set_active_element(e);

			ex_sln.set_active_element(e);

			Quad3D *quad = get_quadrature(e->get_mode());
			ex_sln.set_quad(quad);

			// test values
			if (!(passed &= test_vertex_values(ex_sln, rm, quad))) break;
			for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++)
				if (!(passed &= test_edge_values(iedge, order.get_edge_order(iedge), ex_sln, rm, quad))) break;
			for (int iface = 0; iface < Hex::NUM_FACES; iface++)
				if (!(passed &= test_face_values(iface, order.get_face_order(iface), ex_sln, rm, quad))) break;
			if (!(passed &= test_elem_values(order, ex_sln, rm, quad))) break;
		}

		if (!passed) break;				// do not continue on the error
	}

	(passed) ? printf("Ok\n") : printf("Failed\n");
	if (!passed) res = ERR_FAILURE;

	return res;
}
