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

	forder_t operator+(const forder_t &o) { return max(this->order, o.order); }
	forder_t operator*(const forder_t &o) { return this->order + o.order; }
	// the rest of possible operators

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


///////////////////////////////


// Example of bilinear form definition
template<typename f_t, typename res_t>
res_t int_grad_u_grad_v(int n, double *wt, f_t *u, f_t *v) {
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i]);
	return result;
}

// Example of linear form definition
template<typename f_t, typename res_t>
res_t int_x_dx(int n, double *wt, f_t *u) {
	res_t result = 0.0;
	return result;
}

// Example of surface bilinear form definition
template<typename f_t, typename res_t>
res_t surf_int_u_v(int n, double *wt, f_t *u, f_t *v, int marker) {
	res_t result = 0.0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->fn[i] * v->fn[i] + u->fn[i] * v->fn[i] + u->fn[i] * v->fn[i]);
	return result;
}

// TODO: Example of surface linear form definition


// Testing

void grad_u_grad_v() {
	double fake_wt = 1.0;		// has to be 1.0 (not to ruin mask and order), we do: result += wt[i] * <expression>
	// orders will be taken from PrecalcShapeset class
	fn_order_t order_u(2), order_v(4);
	int order = int_grad_u_grad_v<fn_order_t, int>(1, &fake_wt, &order_u, &order_v);
	printf("Order = %d\n", order);

	// get the mask to know what to precalculate
	fn_mask_t mask_u, mask_v;
	int mask = int_grad_u_grad_v<fn_mask_t, int>(1, &fake_wt, &mask_u, &mask_v);
	printf("Mask =");
	if (mask & FN_VAL) printf(" VAL");
	if (mask & FN_DX) printf(" DX");
	if (mask & FN_DY) printf(" DY");
	if (mask & FN_DZ) printf(" DZ");
	if (mask & FN_DXY) printf(" DXY");
	if (mask & FN_DYZ) printf(" DYZ");
	if (mask & FN_DXZ) printf(" DXZ");
	printf("\n");

	// get the values and transform them to phys. domain (for now just some stupid numbers)
	double wt[] = { 1.0, 2.0, 1.0 };
	double dudx[] = { 0.5, 0.5, 0.5 };
	double dudy[] = { 1.5, 1.5, 1.5 };
	double dudz[] = { 2.5, 2.5, 2.5 };

	fn_t u, v;
	u.dx = dudx; u.dy = dudy; u.dz = dudz;
	v.dx = dudx; v.dy = dudy; v.dz = dudz;
	double result = int_grad_u_grad_v<fn_t, double>(3, wt, &u, &v);
	printf("result = %lf\n", result);
}

int main(int argc, char *argv[]) {
	grad_u_grad_v();

	return -1;
}
