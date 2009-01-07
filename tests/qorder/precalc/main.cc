//
// main.cc
//
// Testing if the values calculated by PrecalcShapset class are correct.
// The purpose is to check if the stuff using Qorder works okey.
// Using H1 shapeset
//
// This is run only with a mesh that is identical to the ref. domain. Otherwise,
// we whould have to transform derivatives (but we test this in qorder-solution)
//
// TODO:
// - test transformations
// - test on vector-valued shapesets
//

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

#define EPS											10e-13


bool test_vertex_values(PrecalcShapeset &pss, Quad3D *&quad) {
	pss.set_quad_order(VTX_QORDER());
	Shapeset *ss = pss.get_shapeset();
	int index = pss.get_active_shape();

	double *val = pss.get_fn_values();
	double *dx = pss.get_dx_values();
	double *dy = pss.get_dy_values();
	double *dz = pss.get_dz_values();
	int np = quad->get_vertex_num_points();
	QuadPt3D *pt = quad->get_vertex_points();

	for (int k = 0; k < np; k++) {
		if (fabs(val[k] - ss->get_fn_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;

		if (fabs(dx[k] - ss->get_dx_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
		if (fabs(dy[k] - ss->get_dy_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
		if (fabs(dz[k] - ss->get_dz_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
	}

	return true;
}

bool test_edge_values(int edge, order1_t order, PrecalcShapeset &pss, Quad3D *&quad) {
	pss.set_quad_order(EDGE_QORDER(edge, order));
	Shapeset *ss = pss.get_shapeset();
	int index = pss.get_active_shape();

	double *val = pss.get_fn_values();
	double *dx = pss.get_dx_values();
	double *dy = pss.get_dy_values();
	double *dz = pss.get_dz_values();
	int np = quad->get_edge_num_points(order);
	QuadPt3D *pt = quad->get_edge_points(edge, order);

	for (int k = 0; k < np; k++) {
		if (fabs(val[k] - ss->get_fn_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;

		if (fabs(dx[k] - ss->get_dx_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
		if (fabs(dy[k] - ss->get_dy_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
		if (fabs(dz[k] - ss->get_dz_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
	}

	return true;
}

bool test_face_values(int face, order2_t order, PrecalcShapeset &pss, Quad3D *&quad) {
	pss.set_quad_order(FACE_QORDER(face, order));
	Shapeset *ss = pss.get_shapeset();
	int index = pss.get_active_shape();

	double *val = pss.get_fn_values();
	double *dx = pss.get_dx_values();
	double *dy = pss.get_dy_values();
	double *dz = pss.get_dz_values();
	int np = quad->get_face_num_points(face, order);
	QuadPt3D *pt = quad->get_face_points(face, order);

	for (int k = 0; k < np; k++) {
		if (fabs(val[k] - ss->get_fn_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;

		if (fabs(dx[k] - ss->get_dx_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
		if (fabs(dy[k] - ss->get_dy_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
		if (fabs(dz[k] - ss->get_dz_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
	}

	return true;
}

bool test_elem_values(order3_t order, PrecalcShapeset &pss, Quad3D *&quad) {
	pss.set_quad_order(ELEM_QORDER(order));
	Shapeset *ss = pss.get_shapeset();
	int index = pss.get_active_shape();

	double *val = pss.get_fn_values();
	double *dx = pss.get_dx_values();
	double *dy = pss.get_dy_values();
	double *dz = pss.get_dz_values();
	int np = quad->get_num_points(order);
	QuadPt3D *pt = quad->get_points(order);

	for (int k = 0; k < np; k++) {
		if (fabs(val[k] - ss->get_fn_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;

		if (fabs(dx[k] - ss->get_dx_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
		if (fabs(dy[k] - ss->get_dy_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
		if (fabs(dz[k] - ss->get_dz_value(index, pt[k].x, pt[k].y, pt[k].z, 0)) > EPS) return false;
	}

	return true;
}

bool test_values(order3_t order, PrecalcShapeset &pss, Quad3D *&quad) {
	// test values
	if (!test_vertex_values(pss, quad)) return false;
	for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++)
		if (!test_edge_values(iedge, order.get_edge_order(iedge), pss, quad)) return false;
	for (int iface = 0; iface < Hex::NUM_FACES; iface++)
		if (!test_face_values(iface, order.get_face_order(iface), pss, quad)) return false;
	if (!test_elem_values(order, pss, quad)) return false;

	return true;
}

// Main ////

int main(int argc, char *argv[]) {
	int res = ERROR_SUCCESS;

	if (argc < 2) return ERR_NOT_ENOUGH_PARAMS;

	H1ShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", argv[1]);
		return ERR_FAILURE;
	}

	// start testing
	bool passed = true;
	for (int d = 1; d <= MAX_QUAD_ORDER; d++) {
		order3_t order(d, d, d);			// quad order
		printf("Order = %d\n", d);

		FOR_ALL_ACTIVE_ELEMENTS(idx, &mesh) {
			Element *e = mesh.elements[idx];

			Quad3D *quad = get_quadrature(e->get_mode());

			pss.set_quad(quad);
			pss.set_active_element(e);

			// vertex functions
			for (int i = 0; i < Hex::NUM_VERTICES; i++) {
				int index = shapeset.get_vertex_index(i);
				pss.set_active_shape(index);
				if (!(passed &= test_values(order, pss, quad))) goto exit;
			}
			// edge functions
			for (int iedge = 0; iedge < Hex::NUM_EDGES; iedge++) {
				int eorder = MAX_ELEMENT_ORDER;
				int *indices = shapeset.get_edge_indices(iedge, 0, eorder);
				for (int i = 0; i < shapeset.get_num_edge_fns(eorder); i++) {
					int index = indices[i];
					pss.set_active_shape(index);
					if (!(passed &= test_values(order, pss, quad))) goto exit;
				}
			}
			// face functions
			for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
				order2_t forder(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
				int *indices = shapeset.get_face_indices(iface, 0, forder);
				for (int i = 0; i < shapeset.get_num_face_fns(forder); i++) {
					int index = indices[i];
					pss.set_active_shape(index);
					if (!(passed &= test_values(order, pss, quad))) goto exit;
				}
			}
			// bubble functions
			order3_t o(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
			int *indices = shapeset.get_bubble_indices(o);
			for (int i = 0; i < shapeset.get_num_bubble_fns(o); i++) {
				int index = indices[i];
				pss.set_active_shape(index);
				if (!(passed &= test_values(order, pss, quad))) goto exit;
			}
		}
	}

exit:
	(passed) ? printf("Ok\n") : printf("Failed\n");
	if (!passed) res == ERR_FAILURE;

	return res;
}
