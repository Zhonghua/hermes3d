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

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

inline scalar integral_u(int o, RealFunction *fu, RefMap *ru) {
	Quad3D *quad = fu->get_quad();

	assert(quad->get_mode() == MODE_HEXAHEDRON);
//	// limit of the order is handled in add_hex_orders function
//	int o = add_hex_orders(fu->get_fn_order(), ru->get_inv_ref_order());
//	LIMIT_HEX_ORDER(o);

	fu->set_quad_order(o);
	double *uval = fu->get_fn_values();
	H1_INTEGRATE_EXPRESSION(uval[i]);
	return result;
}

/// Integral \u \v
///
/// @ingroup h1intergrals
inline scalar integral_u_v(int o, RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();

//	int o = 0;
	assert(quad->get_mode() == MODE_HEXAHEDRON);
//	// limit of the order is handled in add_hex_orders function
//	o = add_hex_orders(fu->get_fn_order(), fv->get_fn_order());
//	o = add_hex_orders(o, ru->get_inv_ref_order());
//	LIMIT_HEX_ORDER(o);

	fu->set_quad_order(o);
	fv->set_quad_order(o);

	double *uval = fu->get_fn_values();
	double *vval = fv->get_fn_values();

	H1_INTEGRATE_EXPRESSION(uval[i] * vval[i]);
	return result;
}

inline scalar integral_grad_u_grad_v(int o, RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();

	assert(quad->get_mode() == MODE_HEXAHEDRON);
/*	int o = 0;

	// limit of the order is handled in add_hex_orders function
	o = add_hex_orders(fu->get_fn_order(), fv->get_fn_order());
	o = add_hex_orders(o, ru->get_inv_ref_order());
//	o = add_hex_orders(o, ru->get_inv_ref_order());
//	o = MAKE_HEX_ORDER(MAX_QUAD_ORDER, MAX_QUAD_ORDER, MAX_QUAD_ORDER);
	LIMIT_HEX_ORDER(o);
//	printf("- : (%d, %d, %d) + (%d, %d, %d) + (%d, %d, %d) => (%d, %d, %d)\n",
//	    GET_HEX_ORDER_1(fu->get_fn_order()), GET_HEX_ORDER_2(fu->get_fn_order()), GET_HEX_ORDER_3(fu->get_fn_order()),
//	    GET_HEX_ORDER_1(fv->get_fn_order()), GET_HEX_ORDER_2(fv->get_fn_order()), GET_HEX_ORDER_3(fv->get_fn_order()),
//	    GET_HEX_ORDER_1(ru->get_inv_ref_order()), GET_HEX_ORDER_2(ru->get_inv_ref_order()), GET_HEX_ORDER_3(ru->get_inv_ref_order()),
//	    GET_HEX_ORDER_1(o), GET_HEX_ORDER_2(o), GET_HEX_ORDER_3(o));
*/
	fu->set_quad_order(o);
	fv->set_quad_order(o);

	double *dudx, *dudy, *dudz;
	fu->get_dx_dy_dz_values(dudx, dudy, dudz);
	double *dvdx, *dvdy, *dvdz;
	fv->get_dx_dy_dz_values(dvdx, dvdy, dvdz);

	H1_INTEGRATE_DD_EXPRESSION(T_DUDX * T_DVDX + T_DUDY * T_DVDY + T_DUDZ * T_DVDZ);
	return result;
}

// surface integrals

inline scalar surf_integral_v(int order, RealFunction *fv, RefMap *rv, FacePos *fp) {
	Quad3D *quad = fv->get_quad();

	assert(quad->get_mode() == MODE_HEXAHEDRON);
	int face_order = 0;
	// limit of the order is handled in add_hex_orders function
//	face_order = add_hex_orders(fv->get_fn_order(), rv->get_inv_ref_order());
//	LIMIT_HEX_ORDER(face_order);
	face_order = get_hex_face_order(fp->face, order);
	int o = MAKE_FACE_ORDER(fp->face, face_order);

	fv->set_quad_order(o, FN_VAL);

	double *vval = fv->get_fn_values();

	H1_INTEGRATE_SURF_EXPRESSION(vval[i]);
	return result;
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

	Element *e = mesh.elements[0];
	pss.set_quad(get_quadrature(e->get_mode()));
	spss.set_quad(get_quadrature(e->get_mode()));

	spss.set_active_element(e);
//	spss->set_master_transform();

    refmap.set_active_element(e);
	refmap.force_transform(pss.get_transform(), pss.get_ctm());

	PrecalcShapeset *fv = &spss;
	PrecalcShapeset *fu = &pss;

//	int iidx = shapeset.get_vertex_index(7);
	const int *edge0_idx = shapeset.get_edge_indices(0, 0, 10);
	const int *edge1_idx = shapeset.get_edge_indices(0, 1, 10);

	const int *face0_idx = shapeset.get_face_indices(2, 0, MAKE_QUAD_ORDER(10, 10));
	const int *face1_idx = shapeset.get_face_indices(2, 4, MAKE_QUAD_ORDER(10, 10));

//	int iidx = shapeset.get_vertex_index(7);
	int iidx = face1_idx[0];
//	int jidx = shapeset.get_vertex_index(7);
	int jidx = face0_idx[0];
	fu->set_active_shape(jidx);
	fv->set_active_shape(iidx);

	RefMap *ru = &refmap;
	RefMap *rv = &refmap;

	int uorder = shapeset.get_order(jidx);
	int vorder = shapeset.get_order(iidx);

//	printf("order: fu = (%d, %d, %d), fv = (%d, %d, %d)\n",
//		GET_HEX_ORDER_1(uorder), GET_HEX_ORDER_2(uorder), GET_HEX_ORDER_3(uorder),
//		GET_HEX_ORDER_1(vorder), GET_HEX_ORDER_2(vorder), GET_HEX_ORDER_3(vorder));

	printf("order: fu = (%d, %d, %d)\n", GET_HEX_ORDER_1(uorder), GET_HEX_ORDER_2(uorder), GET_HEX_ORDER_3(uorder));
	printf("order: fv = (%d, %d, %d)\n", GET_HEX_ORDER_1(vorder), GET_HEX_ORDER_2(vorder), GET_HEX_ORDER_3(vorder));

	int order;
	double val, v;
	double ref;


#ifndef VOL
	order = MAKE_HEX_ORDER(23, 23, 23);
	ref = val = integral_u(order, fu, ru);
	v = int_u(fu, ru);
	printf("order (%d, %d, %d): val = %lf | v = %lf\n",
		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val, v);

	order = uorder;
	val = integral_u(order, fu, ru);
	v = int_u(fu, ru);
	printf("order (%d, %d, %d): val = %lf | err = %e\n",
		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val, fabs(ref - val));

	for (int i = 1; i < 24; i++) {
		order = MAKE_HEX_ORDER(2, 1, i);
		val = integral_u(order, fu, ru);
		printf("order (%d, %d, %d): val = %lf | err = %e\n",
			GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val, fabs(ref - val));
	}

	//
	order = add_hex_orders(uorder, vorder);
	ref = val = integral_u_v(order, fu, fv, ru, rv);
	v = int_u_v(fu, fv, ru, rv);
	printf("order (%d, %d, %d): val = %lf | v = %lf\n",
		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val, v);

	for (int i = 1; i < 24; i++) {
		order = MAKE_HEX_ORDER(4, 2, i);
		val = integral_u_v(order, fu, fv, ru, rv);
		v = int_u_v(fu, fv, ru, rv);
		printf("order (%d, %d, %d): val = %lf | v = %lf | %e\n",
			GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val, v, fabs(ref - val));
	}

	order = MAKE_HEX_ORDER(24, 24, 24);
	val = integral_u_v(order, fu, fv, ru, rv);
	printf("order (%d, %d, %d): val = %lf\n",
		GET_HEX_ORDER_1(order), GET_HEX_ORDER_2(order), GET_HEX_ORDER_3(order), val);

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

	order = MAKE_HEX_ORDER(24, 24, 24);
	printf("%d, %X\n", order, order);
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
