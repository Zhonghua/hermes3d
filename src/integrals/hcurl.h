#ifndef _INTEGRALS_HCURL_H_
#define _INTEGRALS_HCURL_H_

#include "refmap.h"
#include <common/trace.h>
#include <common/error.h>

/// @defgroup hcurlintergrals HCurl intergals

#define HCURL_INTEGRATE_EXPRESSION(exp) \
	scalar result = 0.0; \
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
		mu = ru->get_inv_ref_map(qord); \
		mv = rv->get_inv_ref_map(qord); \
		double *jac = ru->get_jacobian(qord); \
		for (int i = 0; i < np; i++, mu++, mv++) \
			result += pt[i].w * jac[i] * (exp); \
	}}

#define HCURL_INTEGRATE_CURL_EXPRESSION(exp) \
	scalar result = 0.0; \
	{ QuadPt3D *pt = quad->get_points(o); \
	int np = quad->get_num_points(o); \
	double3x3 *mv, *mu; \
	if (ru->is_jacobian_const()) { \
		mu = ru->get_const_ref_map(); \
		mv = rv->get_const_ref_map(); \
		for (int i = 0; i < np; i++) \
			result += pt[i].w * (exp); \
		result /= ru->get_const_jacobian(); \
	} \
	else { \
		mu = ru->get_ref_map(qord); \
		mv = rv->get_ref_map(qord); \
		double *jac = ru->get_jacobian(qord); \
		for (int i = 0; i < np; i++, mu++, mv++) \
			result += pt[i].w  * (exp) / jac[i]; \
	}}

#define HCURL_INTEGRATE_SURF_EXPRESSION(exp) \
	scalar result = 0.0; \
	{ int np = quad->get_face_num_points(fp->face, face_order); \
	fp->space = fp->space_v; \
	QuadPt3D *pt = quad->get_face_points(fp->face, face_order); \
	Point3D *nu, *nv; \
	double3x3 *mv, *mu; \
	if (rv->is_face_const_jacobian(fp->face)) { \
		nu = ru->get_face_const_normal(fp->face); \
		nv = rv->get_face_const_normal(fp->face); \
		mu = ru->get_face_const_inv_ref_map(fp->face); \
		mv = rv->get_face_const_inv_ref_map(fp->face); \
		for (int i = 0; i < np; i++) \
			result += pt[i].w * (exp); \
		result *= rv->get_face_const_jacobian(fp->face); \
	} \
	else { \
		double *jac = rv->get_face_jacobian(fp->face, face_order); \
		nu = ru->get_face_normal(fp->face, face_order); \
		nv = rv->get_face_normal(fp->face, face_order); \
		mu = ru->get_face_inv_ref_map(fp->face, face_order); \
		mv = rv->get_face_inv_ref_map(fp->face, face_order); \
		for (int i = 0; i < np; i++, nu++, nv++, mu++, mv++){ \
			result += pt[i].w * jac[i] * (exp); \
			/*printf("integ point %d, jac %lf, normal (%lf, %lf, %lf), exp " SPS " \n", i, jac[i], (*nu).x, (*nu).y, (*nu).z, SP(std::complex<double>(exp)));*/\
		}\
}}


#define U_CURL_0 (du2dy[i] - du1dz[i])
#define U_CURL_1 (du0dz[i] - du2dx[i])
#define U_CURL_2 (du1dx[i] - du0dy[i])

#define V_CURL_0 (dv2dy[i] - dv1dz[i])
#define V_CURL_1 (dv0dz[i] - dv2dx[i])
#define V_CURL_2 (dv1dx[i] - dv0dy[i])

#define T_U_CURL_0 (U_CURL_0 * (*mu)[0][0] + U_CURL_1 * (*mu)[0][1] + U_CURL_2 * (*mu)[0][2])
#define T_U_CURL_1 (U_CURL_0 * (*mu)[1][0] + U_CURL_1 * (*mu)[1][1] + U_CURL_2 * (*mu)[1][2])
#define T_U_CURL_2 (U_CURL_0 * (*mu)[2][0] + U_CURL_1 * (*mu)[2][1] + U_CURL_2 * (*mu)[2][2])

#define T_V_CURL_0 (V_CURL_0 * (*mv)[0][0] + V_CURL_1 * (*mv)[0][1] + V_CURL_2 * (*mv)[0][2])
#define T_V_CURL_1 (V_CURL_0 * (*mv)[1][0] + V_CURL_1 * (*mv)[1][1] + V_CURL_2 * (*mv)[1][2])
#define T_V_CURL_2 (V_CURL_0 * (*mv)[2][0] + V_CURL_1 * (*mv)[2][1] + V_CURL_2 * (*mv)[2][2])

#define T_U_0 (u0[i] * (*mu)[0][0] + u1[i] * (*mu)[0][1] + u2[i] * (*mu)[0][2])
#define T_U_1 (u0[i] * (*mu)[1][0] + u1[i] * (*mu)[1][1] + u2[i] * (*mu)[1][2])
#define T_U_2 (u0[i] * (*mu)[2][0] + u1[i] * (*mu)[2][1] + u2[i] * (*mu)[2][2])

#define T_V_0 (v0[i] * (*mv)[0][0] + v1[i] * (*mv)[0][1] + v2[i] * (*mv)[0][2])
#define T_V_1 (v0[i] * (*mv)[1][0] + v1[i] * (*mv)[1][1] + v2[i] * (*mv)[1][2])
#define T_V_2 (v0[i] * (*mv)[2][0] + v1[i] * (*mv)[2][1] + v2[i] * (*mv)[2][2])

//for surface integrals : unit normal x u
#define N_X_U_0 ((*nu).y * T_U_2 - (*nu).z * T_U_1)
#define N_X_U_1 ((*nu).z * T_U_0 - (*nu).x * T_U_2)
#define N_X_U_2 ((*nu).x * T_U_1 - (*nu).y * T_U_0)

//for surface integrals : unit normal x v
#define N_X_V_0 ((*nv).y * T_V_2 - (*nv).z * T_V_1)
#define N_X_V_1 ((*nv).z * T_V_0 - (*nv).x * T_V_2)
#define N_X_V_2 ((*nv).x * T_V_1 - (*nv).y * T_V_0)

//for surface integrals : unit normal x u x unit normal ... tangent proejction of u
#define N_X_U_X_N_0 (N_X_U_1 * (*nu).z - N_X_U_2 * (*nu).y)
#define N_X_U_X_N_1 (N_X_U_2 * (*nu).x - N_X_U_0 * (*nu).z)
#define N_X_U_X_N_2 (N_X_U_0 * (*nu).y - N_X_U_1 * (*nu).x)

//for surface integrals : unit normal x v x unit normal ... tangent proejction of v
#define N_X_V_X_N_0 (N_X_V_1 * (*nv).z - N_X_V_2 * (*nv).y)
#define N_X_V_X_N_1 (N_X_V_2 * (*nv).x - N_X_V_0 * (*nv).z)
#define N_X_V_X_N_2 (N_X_V_0 * (*nv).y - N_X_V_1 * (*nv).x)

/// Integral \u \v
///
/// @ingroup hcurlintergrals
inline scalar hcurl_int_u_v(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();

	order3_t o = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	fv->set_quad_order(qord);

	double *u0 = fu->get_fn_values(0);
	double *u1 = fu->get_fn_values(1);
	double *u2 = fu->get_fn_values(2);
	double *v0 = fv->get_fn_values(0);
	double *v1 = fv->get_fn_values(1);
	double *v2 = fv->get_fn_values(2);

	HCURL_INTEGRATE_EXPRESSION(T_U_0 * T_V_0 + T_U_1 * T_V_1 + T_U_2 * T_V_2);
	return result;
}

/// Integral \curl u \curl v
///
/// @ingroup hcurlintergrals
inline scalar hcurl_int_curl_u_curl_v(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv) {
	Quad3D *quad = fu->get_quad();

	order3_t o = fu->get_fn_order() + fv->get_fn_order() + ru->get_ref_order();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);
	fv->set_quad_order(qord);

	double *du0dx, *du0dy, *du0dz, *du1dx, *du1dy, *du1dz, *du2dx, *du2dy, *du2dz;
	double *dv0dx, *dv0dy, *dv0dz, *dv1dx, *dv1dy, *dv1dz, *dv2dx, *dv2dy, *dv2dz;

	fu->get_dx_dy_dz_values(du0dx, du0dy, du0dz, 0);
	fu->get_dx_dy_dz_values(du1dx, du1dy, du1dz, 1);
	fu->get_dx_dy_dz_values(du2dx, du2dy, du2dz, 2);

	fv->get_dx_dy_dz_values(dv0dx, dv0dy, dv0dz, 0);
	fv->get_dx_dy_dz_values(dv1dx, dv1dy, dv1dz, 1);
	fv->get_dx_dy_dz_values(dv2dx, dv2dy, dv2dz, 2);

	HCURL_INTEGRATE_CURL_EXPRESSION(T_U_CURL_0 * T_V_CURL_0 + T_U_CURL_1 * T_V_CURL_1 + T_U_CURL_2 * T_V_CURL_2);
	return result;
}

/// Integral \F \u
///
/// @ingroup hcurlintergrals
inline scalar hcurl_int_F_u(scalar (*F)(double x, double y, double z, int comp), RealFunction *fu, RefMap *ru) {
	Quad3D *quad = fu->get_quad();

	order3_t o(0, 0, 0);
	o.set_maximal();
	qorder_t qord = ELEM_QORDER(o);
	fu->set_quad_order(qord);

	double *u0 = fu->get_fn_values(0);
	double *u1 = fu->get_fn_values(1);
	double *u2 = fu->get_fn_values(2);

	double *x = ru->get_phys_x(o);
	double *y = ru->get_phys_y(o);
	double *z = ru->get_phys_z(o);

	//just to make macro work
	RefMap *rv = ru;

	HCURL_INTEGRATE_EXPRESSION(T_U_0 * F(x[i], y[i], z[i], 0) + T_U_1 * F(x[i], y[i], z[i], 1) + T_U_2 * F(x[i], y[i], z[i], 2));
	return result;
}

#define hcurl_int_F_v(F, fv, rv) hcurl_int_F_u(F, fv, rv)


// surface integrals //////////////////////////////////////////////////////////////////////////////

/// Integral \G dot \v projection (surface)
///
/// @ingroup hcurlintergrals
inline scalar hcurl_surf_int_G_v(RealFunction *fv, RefMap *rv, FacePos *fp) {
	Quad3D *quad = fv->get_quad();

	order3_t order(0, 0, 0);
	order.set_maximal();
	order2_t face_order = order.get_face_order(fp->face);
	qorder_t qord = FACE_QORDER(fp->face, face_order);

	//just to make macro work
	RefMap *ru = rv;

	fv->set_quad_order(qord, FN_VAL);

	double *v0 = fv->get_fn_values(0);
	double *v1 = fv->get_fn_values(1);
	double *v2 = fv->get_fn_values(2);

	double *x = rv->get_face_phys_x(fp->face, face_order);
	double *y = rv->get_face_phys_y(fp->face, face_order);
	double *z = rv->get_face_phys_z(fp->face, face_order);

	HCURL_INTEGRATE_SURF_EXPRESSION(
	        N_X_V_X_N_0 * fp->space->bc_value_callback_by_coord(fp->marker, x[i], y[i], z[i], 0) +
			N_X_V_X_N_1 * fp->space->bc_value_callback_by_coord(fp->marker, x[i], y[i], z[i], 1) +
			N_X_V_X_N_2 * fp->space->bc_value_callback_by_coord(fp->marker, x[i], y[i], z[i], 2));

	return result;
}

/// Integral \u projection dot \v projection (surface)
///
/// @ingroup h1intergrals
inline scalar hcurl_surf_int_u_v(RealFunction *fu, RealFunction *fv, RefMap *ru, RefMap *rv, FacePos *fp) {
	Quad3D *quad = fu->get_quad();

	order3_t order = fu->get_fn_order() + fv->get_fn_order() + ru->get_inv_ref_order();
	order2_t face_order = order.get_face_order(fp->face);
	qorder_t qord = FACE_QORDER(fp->face, face_order);

	fu->set_quad_order(qord, FN_VAL);
	fv->set_quad_order(qord, FN_VAL);

	double *u0 = fu->get_fn_values(0);
	double *u1 = fu->get_fn_values(1);
	double *u2 = fu->get_fn_values(2);
	double *v0 = fv->get_fn_values(0);
	double *v1 = fv->get_fn_values(1);
	double *v2 = fv->get_fn_values(2);

	HCURL_INTEGRATE_SURF_EXPRESSION(N_X_U_X_N_0 * N_X_V_X_N_0 + N_X_U_X_N_1 * N_X_V_X_N_1 + N_X_U_X_N_2 * N_X_V_X_N_2);
	return result;
}


#endif
