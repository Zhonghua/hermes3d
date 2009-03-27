// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2007 - 2009 Pavel Kus <pavel.kus@gmail.com>
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

#ifndef _INTEGRALS_H1_H_
#define _INTEGRALS_H1_H_

#include "refmap.h"
#include "order.h"
#include <common/trace.h>
#include <common/error.h>

/// @defgroup h1intergrals H1 intergals

#define H1_INTEGRATE_EXPRESSION(exp) \
	double result = 0.0; \
	QuadPt3D *pt = quad->get_points(o); \
	int np = quad->get_num_points(o); \
	if (ru->is_jacobian_const()) { \
		for (int i = 0; i < np; i++) \
			result += pt[i].w * (exp); \
		result *= ru->get_const_jacobian(); \
	} \
	else { \
		double *jac = ru->get_jacobian(o); \
		for (int i = 0; i < np; i++) \
			result += pt[i].w * jac[i] * (exp); \
	}


#define H1_INTEGRATE_DD_EXPRESSION(exp) \
	double result = 0.0; \
	{ QuadPt3D *pt = quad->get_points(o); \
	int np = quad->get_num_points(o); \
	double3x3 *mv, *mu; \
	if (ru->is_jacobian_const()) { \
		mu = ru->get_const_inv_ref_map(); \
		mv = rv->get_const_inv_ref_map(); \
		for (int i = 0; i < np; i++) \
			result += pt[i].w * (exp); \
		result *= ru->get_const_jacobian(); \
	} \
	else { \
		mu = ru->get_inv_ref_map(o); \
		mv = rv->get_inv_ref_map(o); \
		double *jac = ru->get_jacobian(o); \
		for (int i = 0; i < np; i++, mu++, mv++) \
			result += pt[i].w * jac[i] * (exp); \
	}}


#define H1_INTEGRATE_SURF_EXPRESSION(exp) \
	scalar result = 0.0; \
	{ int np = quad->get_face_num_points(fp->face, face_order); \
	fp->space = fp->space_v; \
	QuadPt3D *pt = quad->get_face_points(fp->face, face_order); \
	if (rv->is_face_const_jacobian(fp->face)) { \
		for (int i = 0; i < np; i++) \
			result += pt[i].w * (exp); \
		result *= rv->get_face_const_jacobian(fp->face); \
	} \
	else { \
		double *jac = rv->get_face_jacobian(fp->face, face_order); \
		for (int i = 0; i < np; i++) \
			result += pt[i].w * jac[i] * (exp); \
	}}

// inverse matrix is already transposed!
#define T_DUDX (dudx[i] * (*mu)[0][0] + dudy[i] * (*mu)[0][1] + dudz[i] * (*mu)[0][2])
#define T_DUDY (dudx[i] * (*mu)[1][0] + dudy[i] * (*mu)[1][1] + dudz[i] * (*mu)[1][2])
#define T_DUDZ (dudx[i] * (*mu)[2][0] + dudy[i] * (*mu)[2][1] + dudz[i] * (*mu)[2][2])
#define T_DVDX (dvdx[i] * (*mv)[0][0] + dvdy[i] * (*mv)[0][1] + dvdz[i] * (*mv)[0][2])
#define T_DVDY (dvdx[i] * (*mv)[1][0] + dvdy[i] * (*mv)[1][1] + dvdz[i] * (*mv)[1][2])
#define T_DVDZ (dvdx[i] * (*mv)[2][0] + dvdy[i] * (*mv)[2][1] + dvdz[i] * (*mv)[2][2])

/// Integral \u
///
/// @ingroup h1intergrals
inline scalar int_u(RealFunction *fu, RefMap *ru) {
	Quad3D *quad = fu->get_quad();

	order3_t o = fu->get_fn_order() + ru->get_inv_ref_order();
	o.limit();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);

	double *uval = fu->get_fn_values();

	H1_INTEGRATE_EXPRESSION(uval[i]);
	return result;
}

#define int_v(fv, rv) int_u(fv, rv)


/// Integral \u \v
///
/// @ingroup h1intergrals
inline scalar int_u_v(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();

	order3_t o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	fv->set_quad_order(qord);

	double *uval = fu->get_fn_values();
	double *vval = fv->get_fn_values();

	H1_INTEGRATE_EXPRESSION(uval[i] * vval[i]);
	return result;
}


/// Integral \F \u
///
/// @ingroup h1intergrals
inline scalar int_F_u(double (*F)(double x, double y, double z), RealFunction *fu, RefMap *ru) {
	Quad3D *quad = fu->get_quad();

	order3_t o = fu->get_fn_order();
	o.set_maximal();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);

	double *uval = fu->get_fn_values();
	double *x = ru->get_phys_x(o);
	double *y = ru->get_phys_y(o);
	double *z = ru->get_phys_z(o);

	H1_INTEGRATE_EXPRESSION(uval[i] * F(x[i], y[i], z[i]));
	return result;
}

#define int_F_v(F, fv, rv) int_F_u(F, fv, rv)


/// Integral \grad u \grad v
///
/// @ingroup h1intergrals
inline scalar int_grad_u_grad_v(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();

	order3_t o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
	o.limit();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	fv->set_quad_order(qord);

	double *dudx, *dudy, *dudz;
	fu->get_dx_dy_dz_values(dudx, dudy, dudz);
	double *dvdx, *dvdy, *dvdz;
	fv->get_dx_dy_dz_values(dvdx, dvdy, dvdz);

	H1_INTEGRATE_DD_EXPRESSION(T_DUDX * T_DVDX + T_DUDY * T_DVDY + T_DUDZ * T_DVDZ);
	return result;
}

// surface integrals //////////////////////////////////////////////////////////////////////////////

/// Integral \v (surface)
///
/// @ingroup h1intergrals
inline scalar surf_int_v(RealFunction *fv, RefMap *rv, FacePos *fp) {
	Quad3D *quad = fv->get_quad();

	order3_t order = fv->get_fn_order() + rv->get_inv_ref_order();
	order.limit();
	order2_t face_order = order.get_face_order(fp->face);
	qorder_t qord = FACE_QORDER(fp->face, face_order);
	fv->set_quad_order(qord, FN_VAL);

	double *vval = fv->get_fn_values();

	H1_INTEGRATE_SURF_EXPRESSION(vval[i]);
	return result;
}

/// Integral \G \v (surface)
///
/// @ingroup h1intergrals
inline scalar surf_int_G_v(RealFunction *fv, RefMap *rv, FacePos *fp) {
	Quad3D *quad = fv->get_quad();

	order3_t o = fv->get_fn_order();
	o.set_maximal();
	order2_t face_order = o.get_face_order(fp->face);
	qorder_t qord = FACE_QORDER(fp->face, face_order);
	fv->set_quad_order(qord, FN_VAL);

	double *vval = fv->get_fn_values();
	double *x = rv->get_face_phys_x(fp->face, face_order);
	double *y = rv->get_face_phys_y(fp->face, face_order);
	double *z = rv->get_face_phys_z(fp->face, face_order);

	H1_INTEGRATE_SURF_EXPRESSION(vval[i] * fp->space->bc_value_callback_by_coord(fp->marker, x[i], y[i], z[i], 0));
	return result;
}

/// Integral \u \v (surface)
///
/// @ingroup h1intergrals
inline scalar surf_int_u_v(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp) {
	Quad3D *quad = fu->get_quad();

	order3_t order = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
	order.limit();
	order2_t face_order = order.get_face_order(fp->face);
	qorder_t qord = FACE_QORDER(fp->face, face_order);

	fu->set_quad_order(qord, FN_VAL);
	fv->set_quad_order(qord, FN_VAL);

	double *uval = fu->get_fn_values();
	double *vval = fv->get_fn_values();

	H1_INTEGRATE_SURF_EXPRESSION(uval[i] * vval[i]);
	return result;
}

//// error & norm integrals ////////////////////////////////////////////////////////////////////////

// TODO: merge this and things in norm.cc

template<typename T>
inline double int_h1_error(Function<T> *fu, Function<T> *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();
	Quad3D *quadv = fv->get_quad();
	assert(quad == fv->get_quad());

	order3_t o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	fv->set_quad_order(qord);

	scalar *fnu = fu->get_fn_values();
	scalar *fnv = fv->get_fn_values();

	scalar *dudx, *dudy, *dudz, *dvdx, *dvdy, *dvdz;
	fu->get_dx_dy_dz_values(dudx, dudy, dudz);
	fv->get_dx_dy_dz_values(dvdx, dvdy, dvdz);

	H1_INTEGRATE_EXPRESSION(sqr(fnu[i] - fnv[i]) + sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]) + sqr(dudz[i] - dvdz[i]));
	return result;
}

template<typename T>
inline double int_h1_semi_error(Function<T> *fu, Function<T> *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();
	Quad3D *quadv = fv->get_quad();
	assert(quad == fv->get_quad());

	order3_t o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	fv->set_quad_order(qord);

	scalar *fnu = fu->get_fn_values();
	scalar *fnv = fv->get_fn_values();

	scalar *dudx, *dudy, *dudz, *dvdx, *dvdy, *dvdz;
	fu->get_dx_dy_values(dudx, dudy, dudz);
	fv->get_dx_dy_values(dvdx, dvdy, dvdz);

	H1_INTEGRATE_EXPRESSION(sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]) + sqr(dudz[i] - dvdz[i]));
	return result;
}

template<typename T>
inline double int_h1_norm(Function<T> *fu, RefMap *ru) {
	Quad3D *quad = fu->get_quad();

	order3_t o = fu->get_fn_order() + ru->get_inv_ref_order();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);

	scalar *fnu = fu->get_fn_values();
	scalar *dudx, *dudy, *dudz;
	fu->get_dx_dy_dz_values(dudx, dudy, dudz);

	H1_INTEGRATE_EXPRESSION(sqr(fnu[i]) + sqr(dudx[i]) + sqr(dudy[i]) + sqr(dudz[i]));
	return result;
}

template<typename T>
inline double int_h1_seminorm(Function<T> *fu, RefMap *ru) {
	Quad3D *quad = fu->get_quad();

	order3_t o = fu->get_fn_order() + ru->get_inv_ref_order();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);

	scalar *fnu = fu->get_fn_values();
	scalar *dudx, *dudy, *dudz;
	fu->get_dx_dy_dz_values(dudx, dudy, dudy);

	H1_INTEGRATE_EXPRESSION(sqr(dudx[i]) + sqr(dudy[i]) + sqr(dudz[i]));
	return result;
}

#endif
