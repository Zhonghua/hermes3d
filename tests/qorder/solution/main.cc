//
// main.cc
//
// Testing if the values calculated by Solution class are correct.
// The purpose is to check if the Qorder stuff works okey.
//
//
// TODO:
// - more functions
//

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

// error should be smaller than this epsilon
#define EPS								10e-10F

double (*exact_solution)(double x, double y, double z, double &dx, double &dy, double &dz);

// FN #1 ////
//

int m = 2, n = 2, o = 2;

double fn1_fnc(double x, double y, double z) {
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

double fn1_dfnc(double x, double y, double z) {
	double ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	double ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	double ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);

	return -(ddxx + ddyy + ddzz);
}

// needed for calculation norms and used by visualizator
double fn1_exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	// u(x, y, z) = pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4)
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fn1_fnc(x, y, z);
}

//

EBCType fn1_bc_types(int marker) {
	return BC_ESSENTIAL;
}

double fn1_bc_values(int marker, double x, double y, double z, int comp) {
	return fn1_fnc(x, y, z);
}

scalar fn1_bilinear_form(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	return int_grad_u_grad_v(fu, fv, ru, rv);
}

scalar fn1_linear_form(RealFunction *fv, RefMap *rv) {
	return int_F_v(fn1_dfnc, fv, rv);
}

////////////////////////

bool test_vertex_values(Solution &sln, ExactSolution &ex, Quad3D *&quad) {
	sln.set_quad_order(VTX_QORDER());

	double *val = sln.get_fn_values();
	int np = quad->get_vertex_num_points();

	RefMap *rm = ex.get_refmap();
	double *x = rm->get_vertex_phys_x();
	double *y = rm->get_vertex_phys_y();
	double *z = rm->get_vertex_phys_z();

	for (int k = 0; k < np; k++) {
		double dx, dy, dz;
		double fn = exact_solution(x[k], y[k], z[k], dx, dy, dz);

		if (fabs(fn - val[k]) > EPS) {
			printf("Inconsistent value at a vertex: diff = % e (k = %d)\n", fabs(val[k] - fn), k);
			return false;
		}
	}

	return true;
}

bool test_edge_values(int edge, int order, Solution &sln, ExactSolution &ex, Quad3D *&quad) {
	sln.set_quad_order(EDGE_QORDER(edge, order));

	double *val = sln.get_fn_values();
	int np = quad->get_edge_num_points(order);

	RefMap *rm = ex.get_refmap();
	double *x = rm->get_edge_phys_x(edge, order);
	double *y = rm->get_edge_phys_y(edge, order);
	double *z = rm->get_edge_phys_z(edge, order);

	for (int k = 0; k < np; k++) {
		double dx, dy, dz;
		double fn = exact_solution(x[k], y[k], z[k], dx, dy, dz);

		if (fabs(fn - val[k]) > EPS) {
			printf("Inconsistent value on edge %d: diff = % e (k = %d)\n", edge, fabs(val[k] - fn), k);
			return false;
		}
	}

	return true;
}

bool test_face_values(int face, int order, Solution &sln, ExactSolution &ex, Quad3D *&quad) {
	sln.set_quad_order(FACE_QORDER(face, order));

	double *val = sln.get_fn_values();
	int np = quad->get_face_num_points(face, order);

	RefMap *rm = ex.get_refmap();
	double *x = rm->get_face_phys_x(face, order);
	double *y = rm->get_face_phys_y(face, order);
	double *z = rm->get_face_phys_z(face, order);

	for (int k = 0; k < np; k++) {
		double dx, dy, dz;
		double fn = exact_solution(x[k], y[k], z[k], dx, dy, dz);

		if (fabs(fn - val[k]) > EPS) {
			printf("Inconsistent value on face %d: diff = % e (k = %d)\n", face, fabs(val[k] - fn), k);
			return false;
		}
	}

	return true;
}

bool test_elem_values(int order, Solution &sln, ExactSolution &ex, Quad3D *&quad) {
	sln.set_quad_order(ELEM_QORDER(order));

	double *val = sln.get_fn_values();
	int np = quad->get_num_points(order);

	RefMap *rm = ex.get_refmap();
	double *x = rm->get_phys_x(order);
	double *y = rm->get_phys_y(order);
	double *z = rm->get_phys_z(order);

	for (int k = 0; k < np; k++) {
		double dx, dy, dz;
		double fn = exact_solution(x[k], y[k], z[k], dx, dy, dz);

		if (fabs(fn - val[k]) > EPS) {
			printf("Inconsistent value in the interior: diff = % e (k = %d)\n", fabs(val[k] - fn), k);
			return false;
		}
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

	H1Space space(&mesh, &shapeset);

	// fn1
	exact_solution = fn1_exact_solution;
	space.set_bc_types(fn1_bc_types);
	space.set_bc_values(fn1_bc_values);

	int order = MAKE_HEX_ORDER(4, 4, 4);
	space.set_uniform_order(order);
	int ndofs = space.assign_dofs();

	UMFPackLinearSolver solver;

	Discretization d(&solver);
	d.set_num_equations(1);
	d.set_spaces(1, &space);
	d.set_pss(1, &pss);

	// fn1
	d.set_bilinear_form(0, 0, fn1_bilinear_form);
	d.set_linear_form(0, fn1_linear_form);

	// assemble siffness matrix
	d.create_stiffness_matrix();
	d.assemble_stiffness_matrix_and_rhs();

	// solve the stiffness matrix
	Solution sln(&mesh);
	bool solved = d.solve_system(1, &sln);


	// test the solution against the exact solution
	// we do NOT use the norm function to avaoid possible problems in them
	bool passed = true;
	ExactSolution ex_sln(&mesh, fn1_exact_solution);
	FOR_ALL_ACTIVE_ELEMENTS(idx, &mesh) {
		Element *e = mesh.elements[idx];

		// use independent ref. map
		ex_sln.set_active_element(e);
		sln.set_active_element(e);

		Quad3D *quad = get_quadrature(e->get_mode());
		sln.set_quad(quad);
		ex_sln.set_quad(quad);

		// test values
		if (!(passed &= test_vertex_values(sln, ex_sln, quad))) break;
		for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++)
			if (!(passed &= test_edge_values(iedge, get_hex_edge_order(iedge, order), sln, ex_sln, quad))) break;
		for (int iface = 0; iface < Hex::NUM_FACES; iface++)
			if (!(passed &= test_face_values(iface, get_hex_face_order(iface, order), sln, ex_sln, quad))) break;
		if (!(passed &= test_elem_values(order, sln, ex_sln, quad))) break;
	}

	(passed) ? printf("Ok\n") : printf("Failed\n");
	if (!passed) res == ERR_FAILURE;

	return res;
}
