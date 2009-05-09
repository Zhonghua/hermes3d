/*
 * forms.cc
 *
 *  Created on: May 15, 2009
 *      Author: andrsd
 */

#include "h3dconfig.h"
#include "forms.h"
#include <common/callstack.h>

void init_geom(geom_t<forder_t> &e) {
	_F_

	static forder_t x[] = { forder_t(order3_t(1, 0, 0)) };
	static forder_t y[] = { forder_t(order3_t(0, 1, 0)) };
	static forder_t z[] = { forder_t(order3_t(0, 0, 1)) };

	static forder_t nx[] = { forder_t(order3_t(1, 0, 0)) };
	static forder_t ny[] = { forder_t(order3_t(0, 1, 0)) };
	static forder_t nz[] = { forder_t(order3_t(0, 0, 1)) };

	static forder_t tx[] = { forder_t(order3_t(1, 0, 0)) };
	static forder_t ty[] = { forder_t(order3_t(0, 1, 0)) };
	static forder_t tz[] = { forder_t(order3_t(0, 0, 1)) };

	e.x = x;   e.y = y;   e.z = z;
	e.nx = nx; e.ny = ny; e.nz = nz;
	e.tx = tx; e.ty = ty; e.tz = tz;
}

void init_geom(geom_t<double> &e, RefMap *rm, const order3_t &order) {
	_F_

	e.x = rm->get_phys_x(order);
	e.y = rm->get_phys_y(order);
	e.z = rm->get_phys_z(order);

	// TODO: normals and tangents
}

void init_geom(geom_t<double> &e, RefMap *rm, int iface, const order2_t &order) {
	_F_

	e.x = rm->get_face_phys_x(iface, order);
	e.y = rm->get_face_phys_y(iface, order);
	e.z = rm->get_face_phys_z(iface, order);

	// TODO: normals and tangents
}

void init_jwt(PrecalcShapeset *fv, RefMap *rm, const order3_t &order, int &np, double *&jwt) {
	_F_

	Quad3D *quad = fv->get_quad();
	QuadPt3D *pt = quad->get_points(order);
	np = quad->get_num_points(order);

	jwt = new double [np];
	if (rm->is_jacobian_const()) {
		double jac = rm->get_const_jacobian();
		for (int i = 0; i < np; i++) jwt[i] = jac * pt[i].w;
	}
	else {
		double *jac = rm->get_jacobian(order);
		for (int i = 0; i < np; i++) jwt[i] = jac[i] * pt[i].w;
	}
}

void init_jwt(PrecalcShapeset *fv, RefMap *rm, int iface, const order2_t &order, int &np, double *&jwt) {
	_F_

	Quad3D *quad = fv->get_quad();
	QuadPt3D *pt = quad->get_face_points(iface, order);
	np = quad->get_face_num_points(iface, order);

	jwt = new double [np];
	if (rm->is_jacobian_const()) {
		double jac = rm->get_face_const_jacobian(iface);
		for (int i = 0; i < np; i++) jwt[i] = jac * pt[i].w;
	}
	else {
		double *jac = rm->get_face_jacobian(iface, order);
		for (int i = 0; i < np; i++) jwt[i] = jac[i] * pt[i].w;
	}
}

void init_fn(PrecalcShapeset *fn, RefMap *rm, const order3_t &order, fn_t &u) {
	_F_

	Quad3D *quad = fn->get_quad();
	int np = quad->get_num_points(order);

	qorder_t qord = ELEM_QORDER(order);
	fn->set_quad_order(qord);

	u.fn = new scalar [np];
	u.dx = new scalar [np];
	u.dy = new scalar [np];
	u.dz = new scalar [np];

	scalar *val = fn->get_fn_values();
	scalar *dudx, *dudy, *dudz;
	fn->get_dx_dy_dz_values(dudx, dudy, dudz);
	double3x3 *m = rm->get_inv_ref_map(order);
	for (int i = 0; i < np; i++) {
		u.fn[i] = val[i];
		u.dx[i] = (dudx[i] * (*m)[0][0] + dudy[i] * (*m)[0][1] + dudz[i] * (*m)[0][2]);
		u.dy[i] = (dudx[i] * (*m)[1][0] + dudy[i] * (*m)[1][1] + dudz[i] * (*m)[1][2]);
		u.dz[i] = (dudx[i] * (*m)[2][0] + dudy[i] * (*m)[2][1] + dudz[i] * (*m)[2][2]);
	}
}

void init_fn(PrecalcShapeset *fn, RefMap *rm, int iface, const order2_t &order, fn_t &u) {
	_F_

	Quad3D *quad = fn->get_quad();
	int np = quad->get_face_num_points(iface, order);

	qorder_t qord = FACE_QORDER(iface, order);
	fn->set_quad_order(qord);

	u.fn = new scalar [np];
	u.dx = new scalar [np];
	u.dy = new scalar [np];
	u.dz = new scalar [np];

	scalar *val = fn->get_fn_values();
	scalar *dudx, *dudy, *dudz;
	fn->get_dx_dy_dz_values(dudx, dudy, dudz);
	double3x3 *m = rm->get_face_inv_ref_map(iface, order);
	for (int i = 0; i < np; i++) {
		u.fn[i] = val[i];
		u.dx[i] = (dudx[i] * (*m)[0][0] + dudy[i] * (*m)[0][1] + dudz[i] * (*m)[0][2]);
		u.dy[i] = (dudx[i] * (*m)[1][0] + dudy[i] * (*m)[1][1] + dudz[i] * (*m)[1][2]);
		u.dz[i] = (dudx[i] * (*m)[2][0] + dudy[i] * (*m)[2][1] + dudz[i] * (*m)[2][2]);
	}
}
