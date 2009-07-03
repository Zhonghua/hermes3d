
#ifndef _QUAD_CHEB_H_
#define _QUAD_CHEB_H_

#include "quad.h"

/// QuadChebTetra is a special "quadrature" consisting of product Chebyshev
/// points on the reference brick. It is used for expressing
/// the solution on an element as a linear combination of monomials.
///
class QuadChebTetra : public Quad3D {
public:
	QuadChebTetra();
	~QuadChebTetra();
};


/// QuadChebHex is a special "quadrature" consisting of product Chebyshev
/// points on the reference brick. It is used for expressing
/// the solution on an element as a linear combination of monomials.
///
class QuadChebHex : public Quad3D {
public:
	QuadChebHex();
	~QuadChebHex();

	virtual QuadPt3D *get_points(const order3_t &order) {
		if (!tables.exists(order.get_idx())) calc_table(order);
		return tables[order.get_idx()];
	}

protected:
	void calc_table(const order3_t &order);
};

#endif
