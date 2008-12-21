/*
 * grad.cc
 *
 * testing correctness of a gradients
 *
 */

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

void h1_int_vol(RealFunction *fu, double3 result) {
	Quad3D *quad = fu->get_quad();

	// integrate with maximum order
	int o = quad->get_max_order();

	fu->set_quad_order(ELEM_QORDER(o), FN_DX | FN_DY | FN_DZ);

	double *dx = fu->get_dx_values();
	double *dy = fu->get_dy_values();
	double *dz = fu->get_dz_values();

	// integrating over reference brick -> jacobian is 1.0 (we do not have to bother with refmap)
	result[0] = result[1] = result[2] = 0.0;
	QuadPt3D *pt = quad->get_points(o);
	int np = quad->get_num_points(o);
	for (int i = 0; i < np; i++) {
		result[0] += pt[i].w * dx[i];
		result[1] += pt[i].w * dy[i];
		result[2] += pt[i].w * dz[i];
	}
}

void h1_int_surf(RealFunction *fu, double3 result) {
	Quad3D *quad = fu->get_quad();

	Point3D norm[] = {
		{ -1.0,  0.0,  0.0 },
		{  1.0,  0.0,  0.0 },
		{  0.0, -1.0,  0.0 },
		{  0.0,  1.0,  0.0 },
		{  0.0,  0.0, -1.0 },
		{  0.0,  0.0,  1.0 },
	};

	// integrating over reference brick -> jacobian is 1.0 (we do not have to bother with refmap)
	result[0] = result[1] = result[2] = 0.0;
	for (int face = 0; face < Hex::NUM_FACES; face++) {
		int face_order = quad->get_face_max_order(face);
		// integrate with maximum order
		qorder_t surf_order = FACE_QORDER(face, face_order);

		fu->set_quad_order(surf_order, FN_VAL);
		double *val = fu->get_fn_values();

		QuadPt3D *pt = quad->get_face_points(face, face_order);
		int np = quad->get_face_num_points(face, face_order);
		for (int i = 0; i < np; i++) {
			result[0] += pt[i].w * norm[face].x * val[i];
			result[1] += pt[i].w * norm[face].y * val[i];
			result[2] += pt[i].w * norm[face].z * val[i];
		}
	}
}

bool test_grad(int fn_idx, Shapeset *shapeset) {
	PrecalcShapeset pss_u(shapeset);
	pss_u.set_quad(get_quadrature(MODE));

	pss_u.set_active_shape(fn_idx);

//	printf("  fn #%d\n", fn_idx);

	printf(".");
	fflush(stdout);			// prevent caching of output (to see that it did not freeze)

	double3 vol_val = { 0.0, 0.0, 0.0 };
	double3 surf_val = { 0.0, 0.0, 0.0 };

	h1_int_vol(&pss_u, vol_val);
	h1_int_surf(&pss_u, surf_val);

//	printf("    % lf == % lf | % g, % d\n", vol_val[0], surf_val[0], fabs(vol_val[0] - surf_val[0]), fabs(vol_val[0] - surf_val[0]) < EPS);
//	printf("    % lf == % lf | % g, % d\n", vol_val[1], surf_val[1], fabs(vol_val[1] - surf_val[1]), fabs(vol_val[1] - surf_val[1]) < EPS);
//	printf("    % lf == % lf | % g, % d\n", vol_val[2], surf_val[2], fabs(vol_val[2] - surf_val[2]), fabs(vol_val[2] - surf_val[2]) < EPS);

	if (fabs(vol_val[0] - surf_val[0]) > EPS || fabs(vol_val[1] - surf_val[1]) > EPS || fabs(vol_val[2] - surf_val[2]) > EPS) {
		printf("\n");
		ERROR("Gradient values for fn #%d do not match", fn_idx);
		return false;
	}

	return true;
}

bool test_gradients(Shapeset *shapeset) {
	printf("IV. gradients\n");

	// vertex fns
	printf("* Vertex functions\n");
	for (int i = 0; i < Hex::NUM_VERTICES; i++) {
		int fn_idx = shapeset->get_vertex_index(i);
		if (!test_grad(fn_idx, shapeset))
			return false;
	}

	// edge fns
	printf("\n* Edge functions\n");
	for (int i = 0; i < Hex::NUM_EDGES; i++) {
		int order = MAX_ELEMENT_ORDER;
		for (int ori = 0; ori < RefHex::get_edge_orientations(); ori++) {
			int *edge_idx = shapeset->get_edge_indices(i, ori, order);
			for (int j = 0; j < shapeset->get_num_edge_fns(order); j++) {
				if (!test_grad(edge_idx[j], shapeset))
					return false;
			}
		}
	}

	// face fns
	printf("\n* Face functions\n");
	for (int i = 0; i < Hex::NUM_FACES; i++) {
		int order = MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
		for (int ori = 0; ori < RefHex::get_face_orientations(i); ori++) {
			int *face_idx = shapeset->get_face_indices(i, ori, order);
			for (int j = 0; j < shapeset->get_num_face_fns(order); j++) {
				if (!test_grad(face_idx[j], shapeset))
					return false;
			}
		}
	}

	// bubble
	printf("\n* Bubble functions\n");
	int order = MAKE_HEX_ORDER(MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER, MAX_ELEMENT_ORDER);
	int *bubble_idx = shapeset->get_bubble_indices(order);
	for (int j = 0; j < shapeset->get_num_bubble_fns(order); j++) {
		if (!test_grad(bubble_idx[j], shapeset))
			return false;
	}

	return true;
}
