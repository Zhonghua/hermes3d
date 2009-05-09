/*
 * main.cc
 *
 *  Created on: Apr 18, 2009
 *      Author: andrsd
 */

// Stuff to work on:
//
// 1) surface integrals - how?
// 2) x, y, z coordinates in the forms
//    probably add a parameter like pt_t *pt, which will be a pointer to array of phys points
//    so we can do pt[i].x, pt[i].y, pt[i].z to get them
// 3) Hcurl integrals in 3D (normals, etc.)
//


#include <hermes3d.h>
#include <iostream>
using namespace std;

/*
typedef double scalar_t;

class fmask_t {
public:
	fmask_t() { mask = 0; }
	fmask_t(int m) { mask = m; }

	// all operators do binary OR
	fmask_t operator+(const fmask_t &o) { return fmask_t(this->mask | o.mask); }
	fmask_t operator*(const fmask_t &o) { return fmask_t(this->mask | o.mask); }
	// the rest of all possible operators

	// We need conversion to int, becuase we do: result += w[i] * fmask_t
	operator int() { return this->mask; }

protected:
	int mask;
};

class fn_mask_t {
public:
	fn_mask_t() {
		fn[0] = FN_VAL;
		dx[0] = FN_DX;
		dy[0] = FN_DY;
		dz[0] = FN_DZ;
	}

	fmask_t fn[1];								// function values
	fmask_t dx[1], dy[1], dz[1];				// derivatives
};


class forder_t {
public:
	forder_t() { order = 0; }
	forder_t(int o) { order = o; }

	forder_t operator+(const forder_t &o) { return forder_t(max(this->order, o.order)); }
	forder_t operator*(const forder_t &o) { return forder_t(this->order + o.order); }
	// the rest of possible operators
	forder_t operator+=(const forder_t &o) {
		this->order = max(this->order, o.order);
		return *this;
	}

	// We need conversion to int, becuase we do: result += w[i] * fmask_t
	operator int() { return this->order; }

protected:
	int order;
};

class fn_order_t {
public:
	fn_order_t(int order) {
		fn[0] = order;
		dx[0] = order;
		dy[0] = order;
		dz[0] = order;
	}

	forder_t fn[1];							// function values
	forder_t dx[1], dy[1], dz[1];			// derivatives
};

// "Function"
//
class fn_t {
public:
	double *fn;						// function values
	double *dx, *dy, *dz;			// derivatives
};

// Geometry
class geom_t {
public:
	double *x, *y, *z;				// coordinates
	double *nx, *ny, *nz;			// normals
	double *tx, *ty, *tz;			// tangents
};

// User defined function
template<typename T>
class user_fn_t {
public:
	T operator()(T x, T y, T z) const { return x*x + y*y + z*z; }
};

*/

/*
template<typename T>
class ud_fn_t : public user_fn_t<T> {
public:
	T operator()(T x, T y, T z) const { return x*x + y*y + z*z; }
};
*/

template<typename T>
T F(T x, T y, T z) {
	return x*x + y*y + z*z;
}


///////////////////////////////

/*
// Example of bilinear form definition
template<typename f_t, typename res_t>
res_t int_grad_u_grad_v(int n, double *wt, f_t *u, f_t *v, geom_t *e) {
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i]);
	return result;
}
*/

// Example of linear form definition
template<typename f_t, typename res_t>
res_t int_x_dx(int n, double *wt, f_t *u, geom_t<res_t> *e) {
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (e->x[i] * u->dx[i]);
	return result;
}

/*
// Example of surface bilinear form definition
template<typename f_t, typename res_t>
res_t surf_int_u_v(int n, double *wt, f_t *u, f_t *v, FacePos *fp, geom_t *e) {
	res_t result = 0.0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->fn[i] * v->fn[i]);
	return result;
}

// Example of surf integral with user defined function
template<typename f_t, typename res_t>
res_t surf_int_G_v(int n, double *wt, user_fn_t<res_t> G, f_t *v, FacePos *fp, geom_t *e) {
	res_t result = 0.0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (G(e->x[i], e->y[i], e->z[i]) * v->fn[i]);
	return result;
}

//H1_INTEGRATE_SURF_EXPRESSION(vval[i] * fp->space->bc_value_callback_by_coord(fp->marker, x[i], y[i], z[i]));

// Example of user-defined function
template<typename f_t, typename res_t>
res_t int_F_u(int n, double *wt, f_t *u, user_fn_t<res_t> F, geom_t *e) {
	res_t result = 0.0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->fn[i] * F(e->x[i], e->y[i], e->z[i]));
	return result;
}
*/

// Testing

void grad_u_grad_v() {
	printf("grad u v\n");

	geom_t<forder_t> oe;

	double fake_wt = 1.0;		// has to be 1.0 (not to ruin mask and order), we do: result += wt[i] * <expression>
	// orders will be taken from PrecalcShapeset class
	fn_order_t order_u(2), order_v(4);
	forder_t order = int_grad_u_grad_v<fn_order_t, forder_t>(1, &fake_wt, &order_u, &order_v, &oe);
//	printf("  Order = (%d, %d, %d)\n", order.order.x, order.order.y, order.order.z);

/*	// get the mask to know what to precalculate
	fn_mask_t mask_u, mask_v;
	int mask = int_grad_u_grad_v<fn_mask_t, int>(1, &fake_wt, &mask_u, &mask_v);
	printf("  Mask =");
	if (mask & FN_VAL) printf(" VAL");
	if (mask & FN_DX) printf(" DX");
	if (mask & FN_DY) printf(" DY");
	if (mask & FN_DZ) printf(" DZ");
	if (mask & FN_DXY) printf(" DXY");
	if (mask & FN_DYZ) printf(" DYZ");
	if (mask & FN_DXZ) printf(" DXZ");
	printf("\n");
*/
	geom_t<double> e;

	// get the values and transform them to phys. domain (for now just some stupid numbers)
	double wt[] = { 1.0, 2.0, 1.0 };
	double dudx[] = { 0.5, 0.5, 0.5 };
	double dudy[] = { 1.5, 1.5, 1.5 };
	double dudz[] = { 2.5, 2.5, 2.5 };

	fn_t u, v;
	u.dx = dudx; u.dy = dudy; u.dz = dudz;
	v.dx = dudx; v.dy = dudy; v.dz = dudz;
	double result = int_grad_u_grad_v<fn_t, double>(3, wt, &u, &v, &e);
	printf("  result = %lf\n", result);
	printf("\n");
}

void ud_fn() {
	printf("User-defined function\n");

	double x[] = { -1, 0, 1};
	double y[] = { 1, 1, 1};
	double z[] = { 1, 1, 1};

	geom_t<double> e;
	e.x = x;
	e.y = y;
	e.z = z;

	forder_t gox[] = { forder_t(order3_t(1, 1, 1)) };
	forder_t goy[] = { forder_t(order3_t(1, 1, 1)) };
	forder_t goz[] = { forder_t(order3_t(1, 1, 1)) };

	geom_t<forder_t> oe;
	oe.x = gox;
	oe.y = goy;
	oe.z = goz;


	double fake_wt = 1.0;		// has to be 1.0 (not to ruin mask and order), we do: result += wt[i] * <expression>
	// orders will be taken from PrecalcShapeset class
	fn_order_t order_u(order3_t(2, 2, 2));
//	ud_fn_t<forder_t> uf;

//	forder_t one(order3_t(1, 1, 1));

//	forder_t order = uf(one, one, one);
//	forder_t order = int_F_v<fn_order_t, forder_t>(1, &fake_wt, &uf, &order_u, &oe);
	forder_t order = int_F_v<fn_order_t, forder_t>(1, &fake_wt, F, &order_u, &oe);
//	printf("  Order = (%d, %d, %d)\n", order.order.x, order.order.y, order.order.z);

	double wt[] = { 1.0, 2.0, 1.0 };
	double fn[] = { 0.5, 0.5, 0.5 };

	fn_t u;
	u.fn = fn;
//	ud_fn_t<double> ufr;

//	double result = int_F_v<fn_t, double>(3, wt, &ufr, &u, &e);
	double result = int_F_v<fn_t, double>(3, wt, F, &u, &e);
	printf("  result = %lf\n", result);
	printf("\n");
}


int main(int argc, char *argv[]) {
//	grad_u_grad_v();
	ud_fn();

//	int i = 2;
//	forder_t fo(2);
//	forder_t r = i * fo;
//	printf("r = %d\n", (int) r);

//	forder_t r = fo * i;
//	printf("r = %d\n", (int) r);

	return -1;
}
