/*
 * forms.h
 *
 *  Created on: May 9, 2009
 *      Author: andrsd
 */

#ifndef FORMS_H_
#define FORMS_H_

#include "quad.h"
#include "function.h"
#include "refmap.h"

// Base type for orders of functions
//
// We defined a special arithmetics with this type to be able to analyze forms
// and determine the necessary integration order.  This will not only work for forms,
// but also work for functions.
class forder_t {
public:
	forder_t() { order = order3_t(0, 0, 0); }
	forder_t(scalar d) { order = order3_t(0, 0, 0); }
	forder_t(order3_t o) { order = o; }

	order3_t get_order() const { return order; }

	forder_t operator+(const forder_t &o) { return forder_t(max(this->order, o.order)); }
	forder_t operator-(const forder_t &o) { return forder_t(max(this->order, o.order)); }
	forder_t operator-(double d) { return *this; }
	forder_t operator*(const forder_t &o) { return forder_t(this->order + o.order); }

	forder_t operator+=(const forder_t &o) { this->order = max(this->order, o.order); return *this; }

protected:
	order3_t order;
};

inline forder_t operator*(const scalar &a, const forder_t &b) { return b; }
inline forder_t operator+(const scalar &a, const forder_t &b) { return b; }
inline forder_t operator-(const scalar &a, const forder_t &b) { return b; }
inline forder_t operator-(const forder_t &a) { return a; }

inline forder_t pow(const forder_t &a, const double &b) { return forder_t((int) ceil(fabs(b)) * a.get_order()); }


///////////////////////////////////////////////////////////////////////////////////////////////////

class fn_order_t {
public:
	fn_order_t(order3_t order) {
		fn[0] = order;
		dx[0] = order;
		dy[0] = order;
		dz[0] = order;
	}

	forder_t fn[1];							// function values
	forder_t dx[1], dy[1], dz[1];			// derivatives
};

// Function
class fn_t {
public:
	double *fn;						// function values
	double *dx, *dy, *dz;			// derivatives

	fn_t() {
		fn = NULL;
		dx = dy = dz = NULL;
	}

	virtual ~fn_t() {
		delete [] fn;
		delete [] dx;
		delete [] dy;
		delete [] dz;
	}
};


// Geometry of the element
template<typename T>
class geom_t {
public:
	T *x, *y, *z;				// coordinates
	T *nx, *ny, *nz;			// normals
	T *tx, *ty, *tz;			// tangents

	geom_t() {
		x = y = z = NULL;
		nx = ny = nz = NULL;
		tx = ty = tz = NULL;
	}
};

void init_geom(geom_t<forder_t> &e);
void init_geom(geom_t<double> &e, RefMap *rm, const order3_t &order);
void init_geom(geom_t<double> &e, RefMap *rm, int iface, const order2_t &order);

void init_jwt(PrecalcShapeset *fv, RefMap *rm, const order3_t &order, int &np, double *&jwt);
void init_jwt(PrecalcShapeset *fv, RefMap *rm, int iface, const order2_t &order, int &np, double *&jwt);
void init_fn (PrecalcShapeset *fu, RefMap *rm, const order3_t &order, fn_t &u);
void init_fn (PrecalcShapeset *fu, RefMap *rm, int iface, const order2_t &order, fn_t &u);


#endif /* FORMS_H_ */
