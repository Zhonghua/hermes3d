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

/*
 * h1-hex-real.cc
 *
 *  Created on: Nov 4, 2008
 *      Author: andrsd
 */

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define EPS											10e-12

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

inline scalar integral_u(order3_t o, RealFunction *fu, RefMap *ru) {
	Quad3D *quad = fu->get_quad();
	assert(quad->get_mode() == MODE_HEXAHEDRON);
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	double *uval = fu->get_fn_values();

	double result = 0.0;
	QuadPt3D *pt = quad->get_points(o);
	int np = quad->get_num_points(o);
	double *jac = ru->get_jacobian(qord);
	for (int i = 0; i < np; i++)
		result += pt[i].w * jac[i] * (uval[i]);

	return result;
}

/// Integral \u \v
///
/// @ingroup h1intergrals
inline scalar integral_u_v(order3_t o, RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();
	assert(quad->get_mode() == MODE_HEXAHEDRON);
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	fv->set_quad_order(qord);

	double *uval = fu->get_fn_values();
	double *vval = fv->get_fn_values();

	double result = 0.0;
	QuadPt3D *pt = quad->get_points(o);
	int np = quad->get_num_points(o);
	double *jac = ru->get_jacobian(qord);
	for (int i = 0; i < np; i++)
		result += pt[i].w * jac[i] * (uval[i] * vval[i]);

	return result;
}

inline scalar integral_grad_u_grad_v(order3_t o, RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();

	assert(quad->get_mode() == MODE_HEXAHEDRON);
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	fv->set_quad_order(qord);

	double *dudx, *dudy, *dudz;
	fu->get_dx_dy_dz_values(dudx, dudy, dudz);
	double *dvdx, *dvdy, *dvdz;
	fv->get_dx_dy_dz_values(dvdx, dvdy, dvdz);

	double result = 0.0;
	QuadPt3D *pt = quad->get_points(o);
	int np = quad->get_num_points(o);
	double3x3 *mv, *mu;
	mu = ru->get_inv_ref_map(qord);
	mv = rv->get_inv_ref_map(qord);
	double *jac = ru->get_jacobian(qord);
	for (int i = 0; i < np; i++, mu++, mv++)
		result += pt[i].w * jac[i] * (T_DUDX * T_DVDX + T_DUDY * T_DVDY + T_DUDZ * T_DVDZ);

	return result;
}

// surface integrals

scalar surf_integral_v(order3_t order, RealFunction *fv, RefMap *rv, FacePos *fp) {
	Quad3D *quad = fv->get_quad();

	assert(quad->get_mode() == MODE_HEXAHEDRON);
	order2_t face_order = order.get_face_order(fp->face);
	qorder_t qord = FACE_QORDER(fp->face, face_order);

	fv->set_quad_order(qord, FN_VAL);
	double *vval = fv->get_fn_values();

	scalar result = 0.0;
	int np = quad->get_face_num_points(fp->face, face_order);
	fp->space = fp->space_v;
	QuadPt3D *pt = quad->get_face_points(fp->face, face_order);
	double *jac = rv->get_face_jacobian(fp->face, face_order);
	for (int i = 0; i < np; i++)
		result += pt[i].w * jac[i] * (vval[i]);

	return result;
}

scalar surf_integral_u_v(order3_t order, RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp) {
	Quad3D *quad = fu->get_quad();

	order2_t face_order = order.get_face_order(fp->face);
	qorder_t qord = FACE_QORDER(fp->face, face_order);

	fu->set_quad_order(qord, FN_VAL);
	fv->set_quad_order(qord, FN_VAL);

	double *uval = fu->get_fn_values();
	double *vval = fv->get_fn_values();

	scalar result = 0.0;
	int np = quad->get_face_num_points(fp->face, face_order);
	fp->space = fp->space_v;
	QuadPt3D *pt = quad->get_face_points(fp->face, face_order);
	double *jac = rv->get_face_jacobian(fp->face, face_order);
	for (int i = 0; i < np; i++)
		result += pt[i].w * jac[i] * (uval[i] * vval[i]);

	return result;
}


void test_vertex_fns(Shapeset *ss, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *refmap) {
	printf("\nTesting vertex functions\n");

	int iidx = ss->get_vertex_index(0);
	int jidx = ss->get_vertex_index(6);

	fu->set_active_shape(jidx);
	fv->set_active_shape(iidx);

	RefMap *ru = refmap;
	RefMap *rv = refmap;

	order3_t uorder = ss->get_order(jidx);
	order3_t vorder = ss->get_order(iidx);

	printf("  fu: order = (%d, %d, %d)\n", uorder.x, uorder.y, uorder.z);
	printf("  fv: order = (%d, %d, %d)\n", vorder.x, vorder.y, vorder.z);

	order3_t order;
	double val;
	double ref;

	// int(u)
	printf("  * int_u\n");
	ref = integral_u(order3_t(23, 23, 23), fu, ru);
	val = int_u(fu, ru);
	printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", uorder.x, uorder.y, uorder.z, val, fabs(ref - val));
	if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

	// int(u v)
	printf("  * int_u_v\n");
	ref = integral_u_v(order3_t(23, 23, 23), fu, fv, ru, rv);
	val = int_u_v(fu, fv, ru, rv);
	order = uorder + vorder;
	printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", order.x, order.y, order.z, val, fabs(ref - val));
	if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

	// grad(u)grad(v)
	printf("  * int_grad_u_grad_v\n");
	ref = integral_grad_u_grad_v(order3_t(23, 23, 23), fu, fv, ru, rv);
	val = int_grad_u_grad_v(fu, fv, ru, rv);
	order = uorder + vorder;
	printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", order.x, order.y, order.z, val, fabs(ref - val));
	if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

	// surf_int_v
	printf("  * surf_int_v\n");
	for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
		FacePos fp;
		fp.face = iface;
		ref = surf_integral_v(order3_t(23, 23, 23), fu, ru, &fp);
		val = surf_int_v(fu, ru, &fp);
		order = uorder;
		printf("    iface = %d; order (%d, %d, %d): val = %lf (diff = %e)\n", iface, order.x, order.y, order.z, val, fabs(ref - val));
		if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }
	}

}

void test_edge_vertex_fns(Shapeset *ss, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *refmap) {
	printf("\nTesting edge/vertex functions\n");

	int o = 10;
	const int *edge0_idx = ss->get_edge_indices(0, 0, o);
	int edge_fns = ss->get_num_edge_fns(o);

	int iidx = ss->get_vertex_index(0);

	for (int j = 2; j < edge_fns; j++) {
		int jidx = edge0_idx[j];

		fu->set_active_shape(jidx);
		fv->set_active_shape(iidx);

		RefMap *ru = refmap;
		RefMap *rv = refmap;

		order3_t uorder = ss->get_order(jidx);
		order3_t vorder = ss->get_order(iidx);

		printf("  fu: order = (%d, %d, %d)\n", uorder.x, uorder.y, uorder.z);
		printf("  fv: order = (%d, %d, %d)\n", vorder.x, vorder.y, vorder.z);

		order3_t order;
		double val;
		double ref;

		// int(u)
		printf("  * int_u\n");
		ref = integral_u(order3_t(23, 23, 23), fu, ru);
		val = int_u(fu, ru);
		printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", uorder.x, uorder.y, uorder.z, val, fabs(ref - val));
		if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

		// int(u v)
		printf("  * int_u_v\n");
		ref = integral_u_v(order3_t(23, 23, 23), fu, fv, ru, rv);
		val = int_u_v(fu, fv, ru, rv);
		order = uorder + vorder;
		printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", order.x, order.y, order.z, val, fabs(ref - val));
		if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

		// grad(u)grad(v)
		printf("  * int_grad_u_grad_v\n");
		ref = integral_grad_u_grad_v(order3_t(23, 23, 23), fu, fv, ru, rv);
		val = int_grad_u_grad_v(fu, fv, ru, rv);
		order = uorder + vorder;
		printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", order.x, order.y, order.z, val, fabs(ref - val));
		if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

		// surf_int_v
		printf("  * surf_int_v\n");
		for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
			FacePos fp;
			fp.face = iface;
			ref = surf_integral_v(order3_t(23, 23, 23), fu, ru, &fp);
			val = surf_int_v(fu, ru, &fp);
			order = uorder;
			printf("    iface = %d; order (%d, %d, %d): val = %lf, ref = %lf (diff = %e)\n", iface, order.x, order.y, order.z, val, ref, fabs(ref - val));
			if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }
		}

		printf("  --\n");
	}
}

void test_edge_edge_fns(Shapeset *ss, PrecalcShapeset *fu, PrecalcShapeset *fv, RefMap *refmap) {
	printf("\nTesting edge/edge functions\n");

	int o = 10;
	const int *edge0_idx = ss->get_edge_indices(0, 0, o);
	const int *edge1_idx = ss->get_edge_indices(1, 0, o);
	int edge_fns = ss->get_num_edge_fns(o);

	for (int j = 0; j < edge_fns; j++) {
		int iidx = edge1_idx[j];

		for (int j = 0; j < edge_fns; j++) {
			int jidx = edge0_idx[j];

			fu->set_active_shape(jidx);
			fv->set_active_shape(iidx);

			RefMap *ru = refmap;
			RefMap *rv = refmap;

			order3_t uorder = ss->get_order(jidx);
			order3_t vorder = ss->get_order(iidx);

			printf("  fu: order = (%d, %d, %d)\n", uorder.x, uorder.y, uorder.z);
			printf("  fv: order = (%d, %d, %d)\n", vorder.x, vorder.y, vorder.z);

			order3_t order;
			double val;
			double ref;

			// int(u)
			printf("  * int_u\n");
			ref = integral_u(order3_t(23, 23, 23), fu, ru);
			val = int_u(fu, ru);
			printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", uorder.x, uorder.y, uorder.z, val, fabs(ref - val));
			if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

			// int(u v)
			printf("  * int_u_v\n");
			ref = integral_u_v(order3_t(23, 23, 23), fu, fv, ru, rv);
			val = int_u_v(fu, fv, ru, rv);
			order = uorder + vorder;
			printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", order.x, order.y, order.z, val, fabs(ref - val));
			if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

			// grad(u)grad(v)
			printf("  * int_grad_u_grad_v\n");
			ref = integral_grad_u_grad_v(order3_t(23, 23, 23), fu, fv, ru, rv);
			val = int_grad_u_grad_v(fu, fv, ru, rv);
			order = uorder + vorder;
			printf("    order (%d, %d, %d): val = %lf (diff = %e)\n", order.x, order.y, order.z, val, fabs(ref - val));
			if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }

			// surf_int_v
			printf("  * surf_int_v\n");
			for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
				FacePos fp;
				fp.face = iface;
				ref = surf_integral_v(order3_t(23, 23, 23), fu, ru, &fp);
				val = surf_int_v(fu, ru, &fp);
				order = uorder;
				printf("    iface = %d; order (%d, %d, %d): val = %lf (diff = %e)\n", iface, order.x, order.y, order.z, val, fabs(ref - val));
				if (fabs(ref - val) > EPS) { throw ERROR_FAILURE; }
			}

			printf("  --\n");
		}
	}
}

int main(int argc, char *argv[]) {
	DEBUG_OUTPUT_OFF;
	SET_VERBOSE_LEVEL(0);

	if (argc < 2) {
		ERROR("Not enough parameters");
		return ERR_NOT_ENOUGH_PARAMS;
	}

	int ret = ERROR_SUCCESS;

	H1ShapesetLobattoHex shapeset;
	PrecalcShapeset pss(&shapeset);

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", argv[1]);
		return ERR_FAILURE;
	}

	RefMap refmap(&mesh);
	PrecalcShapeset spss(&pss);

	assert(mesh.elements.count() > 0);
	Element *e = mesh.elements[1];
	pss.set_quad(get_quadrature(e->get_mode()));
	spss.set_quad(get_quadrature(e->get_mode()));

	spss.set_active_element(e);
//	spss->set_master_transform();

    refmap.set_active_element(e);
	refmap.force_transform(pss.get_transform(), pss.get_ctm());

	PrecalcShapeset *fv = &spss;
	PrecalcShapeset *fu = &pss;

//	int iidx = shapeset.get_vertex_index(7);
//	const int *edge0_idx = shapeset.get_edge_indices(0, 0, 10);
//	const int *edge1_idx = shapeset.get_edge_indices(0, 1, 10);

//	const int *face0_idx = shapeset.get_face_indices(2, 0, order2_t(10, 10));
//	const int *face1_idx = shapeset.get_face_indices(2, 4, order2_t(10, 10));

/*
	int o = 10;
	const int *edge0_idx = shapeset.get_edge_indices(4, 0, o);
	int edge_fns = shapeset.get_num_edge_fns(o);

//	for (int j = 0; j < edge_fns; j++)
	{
		int iidx = edge0_idx[0];

		fu->set_active_shape(iidx);
		RefMap *ru = &refmap;
		order3_t uorder = shapeset.get_order(iidx);

		order3_t order;
		order = order3_t(24, 24, 24);
		double ref = integral_u(order, fu, ru);

		order = order3_t(1, 1, 2);
		double val = integral_u(order, fu, ru);
//		v = int_u(fu, ru);
		printf("order (%d, %d, %d): val = %lf, ref = %lf, %e\n", order.x, order.y, order.z, val, ref, fabs(ref - val));
	}
*/

	try {
		test_vertex_fns(&shapeset, fu, fv, &refmap);
		test_edge_vertex_fns(&shapeset, fu, fv, &refmap);
		test_edge_edge_fns(&shapeset, fu, fv, &refmap);

//		test_face_vertex_fns(&shapeset, fu, fv, &refmap);
//		test_face_edge_fns(&shapeset, fu, fv, &refmap);
//		test_face_face_fns(&shapeset, fu, fv, &refmap);
//		test_bubble_vertex_fns(&shapeset, fu, fv, &refmap);
//		test_bubble_edge_fns(&shapeset, fu, fv, &refmap);
//		test_bubble_face_fns(&shapeset, fu, fv, &refmap);
//		test_bubble_bubble_fns(&shapeset, fu, fv, &refmap);

		printf("Passed\n");
	}
	catch (int e) {
		printf("Failed\n");
		ret = e;
	}

#ifdef VOL
	order = order3_t(23, 23, 23);
	ref = val = integral_u(order, fu, ru);
	v = int_u(fu, ru);
	printf("order (%d, %d, %d): val = %lf | v = %lf\n", order.x, order.y, order.z, val, v);

	order = uorder;
	val = integral_u(order, fu, ru);
	v = int_u(fu, ru);
	printf("order (%d, %d, %d): val = %lf | err = %e\n", order.x, order.y, order.z, val, fabs(ref - val));

	for (int i = 1; i < 24; i++) {
		order = order3_t(2, 1, i);
		val = integral_u(order, fu, ru);
		printf("order (%d, %d, %d): val = %lf | err = %e\n", order.x, order.y, order.z, val, fabs(ref - val));
	}

	//
	order = uorder + vorder;
	ref = val = integral_u_v(order, fu, fv, ru, rv);
	v = int_u_v(fu, fv, ru, rv);
	printf("order (%d, %d, %d): val = %lf | v = %lf\n", order.x, order.y, order.z, val, v);

	for (int i = 1; i < 24; i++) {
		order = order3_t(4, 2, i);
		val = integral_u_v(order, fu, fv, ru, rv);
		v = int_u_v(fu, fv, ru, rv);
		printf("order (%d, %d, %d): val = %lf | v = %lf | %e\n", order.x, order.y, order.z, val, v, fabs(ref - val));
	}

	order = order3_t(24, 24, 24);
	val = integral_u_v(order, fu, fv, ru, rv);
	printf("order (%d, %d, %d): val = %lf\n", order.x, order.y, order.z, val);

//	val = int_u_v(fu, fv, ru, rv);
//	printf("order (%d, %d, %d): val = %lf\n",
//		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val);

/*	//
	order = add_hex_orders(uorder, vorder);
	val = integral_grad_u_grad_v(order, fu, fv, ru, rv);
	printf("order (%d, %d, %d): val = %lf\n",
		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val);

	order = MAKE_HEX_ORDER(24, 24, 24);
	val = integral_grad_u_grad_v(order, fu, fv, ru, rv);
	printf("order (%d, %d, %d): val = %lf\n",
		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val);
*/

//	for (int i = 1; i < 24; i++) {
//		for (int j = 1; j < 24; j++) {
//			for (int k = 1; k < 24; k++) {
//				int order = MAKE_HEX_ORDER(i, j, k);
//
//				double val = integral_grad_u_grad_v(order, fu, fv, ru, rv);
//				printf("order (%d, %d, %d): val = %lf\n",
//					GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val);
//			}
//		}
//	}

/*
	for (int i = 1; i < 24; i++) {
//		for (int j = 1; j < 24; j++) {
//			for (int k = 1; k < 24; k++) {
//				int order = MAKE_HEX_ORDER(i, j, k);
				int order = MAKE_HEX_ORDER(i, 1, 1);

				double val = integral_u(order, fu, ru);
				printf("order (%d, %d, %d): val = %lf\n",
					GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val);
//			}
//		}
	}
*/

//	order = order3_t(24, 24, 24);
//	printf("%d, %X\n", order, order);
#endif

#ifdef SURF
	FacePos fp;
	fp.face = 2;

	order = uorder;
	val = surf_integral_v(order, fu, ru, &fp);
	v = surf_int_v(fu, ru, &fp);
	printf("order (%d, %d, %d): val = %lf | v = %lf\n",
		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val, v);

	for (int i = 1; i < 24; i++) {
		order = MAKE_HEX_ORDER(2, 1, i);
		val = surf_integral_v(order, fu, ru, &fp);
		printf("order (%d, %d, %d): val = %lf\n",
			GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val);
	}

	order = MAKE_HEX_ORDER(23, 23, 23);
	val = surf_integral_v(order, fu, ru, &fp);
	v = surf_int_v(fu, ru, &fp);
	printf("order (%d, %d, %d): val = %lf | v = %lf\n",
		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val, v);
#endif

	return 0;
}
