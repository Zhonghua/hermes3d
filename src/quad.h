#ifndef _QUAD_H_
#define _QUAD_H_

#include "config.h"
#include "common.h"
#include <common/error.h>

/// @defgroup quadratures Numerical quadratures
///
/// There is a class called QuadXD (where X = 1, 2, 3) that represents numerical
/// quadratures in X dimensions. Each class holds quadratures points for every order.
/// The order of integration is bounded by max_order. The class provides interface
/// for getting quadratures points, numbers of points and maximal order of integration.
/// Numerical quadratures of higher dimendions (2, 3) also provides integration points
/// for surface integrals. This simplifies evaluating of these integrals on elements
/// (one do NOT have to take care about transformations, it is done by the QuadXXX
/// class).
///
/// There are prepared quadratures for standard domains (line, triangle, quad, hex,
/// tetra and prism) which are described in Pavel's Gray book (see quadstd.h). You
/// can create your own quadratures by deriving a class from QuadXD.
///
///
///
/// NOTES:
/// - Quad is short name for quadratures, but can be mismatched with Quadrilateral
/// - QuadPtXD - quadratures point (find better name)
///


/// Quadrature point in 1D
///
/// @ingroup quadratures
struct QuadPt1D {
	double x;		// x-coordinate
	double w;		// weight

	QuadPt1D() { }		// default c-tor
	QuadPt1D(double x, double w) {
		this->x = x;
		this->w = w;
	}

	double operator[](int idx) const {
		if (idx == 0) return this->x;
		else { ERROR("Index out of bounds"); return 0; }
	}
};


/// Quadrature point in 2D
///
/// @ingroup quadratures
struct QuadPt2D {
	double x, y;		// x and y-coordinate
	double w;			// weight

	QuadPt2D() { }		// default c-tor
	QuadPt2D(double x, double y, double w) {
		this->x = x;
		this->y = y;
		this->w = w;
	}

	double operator[](int idx) const {
		if (idx == 0) return this->x;
		else if (idx == 1) return this->y;
		else { ERROR("Index out of bounds"); return 0; }
	}
};


/// Quadrature point in 2D
///
/// @ingroup quadratures
struct QuadPt3D {
	double x, y, z;		// x, y and z-coordinate
	double w;			// weight

	QuadPt3D() { } 		// default c-tor
	QuadPt3D(double x, double y, double z, double w) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}

	double operator[](int idx) const {
		if (idx == 0) return this->x;
		else if (idx == 1) return this->y;
		else if (idx == 2) return this->z;
		else { ERROR("Index out of bounds"); return 0; }
	}
};


///
/// 1D quadratures
///

/// General numerical quadratures in 1D
///
/// @ingroup quadratures
class Quad1D {
public:
	QuadPt1D *get_points(int order) const { return tables[order]; }
	int get_num_points(int order) const { return np[order]; };
	int get_max_order() const { return max_order; }

	EMode1D get_mode() const { return mode; }

protected:
	/// mode of quadratures (MODE_LINE)
	EMode1D mode;
	/// maximal order of integration
	int max_order;
	/// table with integration points for each order
	/// indexing: [order][point no.]
	QuadPt1D **tables;
	/// number of integration points for each order
	/// indexing: [order]
	int *np;
};


///
/// 2D quadratures
///


/// Numerical quadratures in 2D
///
/// @ingroup quadratures
class Quad2D {
public:
	QuadPt2D *get_points(int order) const { return tables[order]; }
	inline int get_num_points(int order) const { return np[order]; };

	QuadPt2D *get_edge_points(int edge, int order) const { return edge_tables[edge][order]; }
//	double get_edge_jacobian(int edge) const { return edge_jacobian[edge]; }

	int get_max_order() const { return max_order; }

	EMode2D get_mode() const { return mode; }

protected:
	/// mode of quadratures (MODE_TRIANGLE, MODE_QUAD)
	EMode2D mode;
	/// maximal order for integration (interpretation depened on the mode)
	int max_order;
	/// number of integration points
	/// indexing: [order]
	int *np;
	/// tables with integration points
	/// indexing: [order][point no.]
	QuadPt2D **tables;
	/// tables with integration points for edges (?)
	/// indexing: [edge][order][point no.]
	QuadPt2D ***edge_tables;
//	/// indexing: [edge]
//	double *edge_jacobian;
};


//
// 3D quadratures
//


/// Numerical quadratures in 3D
///
/// @ingroup quadratures
class Quad3D {
public:
	virtual QuadPt3D *get_points(int order) { return tables[order]; }
	virtual int get_num_points(int order) { return np[order]; }

	virtual QuadPt3D *get_edge_points(int edge, int order) { return edge_tables[edge][order]; }
	int get_edge_num_points(int order) const { return np_edge[order]; }
//	double get_edge_jacobian(int edge) const { return edge_jacobian[edge]; }

	virtual QuadPt3D *get_face_points(int face, int order) { return face_tables[face][order]; }
	int get_face_num_points(int face, int order) const { return np_face[order]; }
//	double get_face_jacobian(int face) const { return face_jacobian[face]; }

	virtual QuadPt3D *get_vertex_points() { return vertex_table; }
	int get_vertex_num_points() const { return np_vertex; }

	int get_max_order() const { return max_order; }
	int get_edge_max_order(int edge) const { return max_edge_order; }
	int get_face_max_order(int face) const { return max_face_order; }

	EMode3D get_mode() const { return mode; }

protected:
	/// mode of quadratures (MODE_TETRAHEDRON, MODE_HEXAHEDRON, MODE_PRISM)
	EMode3D mode;
	/// maximal order for integration (interpretation depened on the mode)
	int max_order;
	int max_edge_order;
	int max_face_order;

	QuadPt3D **tables;
	QuadPt3D ***edge_tables;
	QuadPt3D ***face_tables;
	QuadPt3D *vertex_table;
	int *np;
	int *np_edge;
	int *np_face;
	int np_vertex;
//	double *edge_jacobian;
//	double *face_jacobian;
};


// interface for getting quadratures - library wide ////////////////////////////////////////////////

Quad1D *get_quadrature(EMode1D mode);
Quad2D *get_quadrature(EMode2D mode);
Quad3D *get_quadrature(EMode3D mode);


// Helpers for calculating with orders ////////////////////////////////////////////////////////////


// maximal order of quadratures for 1D
#define MAX_QUAD_ORDER								24

// maximal order of quadratures for triangle
#define MAX_QUAD_ORDER_TRI							20
// maximal order of quadratures for tetra
#define MAX_QUAD_ORDER_TETRA						20

// QUAD specific //////////////////////////////////////////////////////////////////////////////////

const int base_quad_coding = MAX_QUAD_ORDER + 1;

#define MAKE_QUAD_ORDER(h_order, v_order) ((v_order) * (base_quad_coding) + (h_order))
#define GET_QUAD_ORDER_1(order) ((order) % (base_quad_coding))
#define GET_QUAD_ORDER_2(order) ((order) / (base_quad_coding))

inline Order2 add_quad_orders(Order2 ord1, Order2 ord2) {
	int o1 = GET_QUAD_ORDER_1(ord1) + GET_QUAD_ORDER_1(ord2);
	if (o1 > MAX_QUAD_ORDER) o1 = MAX_QUAD_ORDER;
	int o2 = GET_QUAD_ORDER_2(ord1) + GET_QUAD_ORDER_2(ord2);
	if (o2 > MAX_QUAD_ORDER) o2 = MAX_QUAD_ORDER;
	return MAKE_QUAD_ORDER(o1, o2);
}

inline void make_quad_debug(Order2 &order) {
	Order1 o1 = GET_QUAD_ORDER_1(order);
	Order1 o2 = GET_QUAD_ORDER_2(order);
	Order1 max = o1 > o2 ? o1 : o2;
	order = MAKE_QUAD_ORDER(max, max);
}

// TETRA specific /////////////////////////////////////////////////////////////////////////////////

#ifndef DEBUG_ORDER
	#define LIMIT_TETRA_ORDER(o) 						if ((o) > MAX_QUAD_ORDER_TETRA) o = MAX_QUAD_ORDER_TETRA;
#else
	#define LIMIT_TETRA_ORDER(o) 						o = MAX_QUAD_ORDER_TETRA;
#endif


// HEX specific ///////////////////////////////////////////////////////////////////////////////////

const int base_hex_coding = MAX_QUAD_ORDER + 1;

#define MAKE_HEX_ORDER(h_order, f_order, v_order) (((v_order) * (base_hex_coding) + (f_order)) * (base_hex_coding) + (h_order))
#define GET_HEX_ORDER_1(order) ((order) % (base_hex_coding))
#define GET_HEX_ORDER_2(order) (((order) / (base_hex_coding)) % (base_hex_coding))
#define GET_HEX_ORDER_3(order) ((order) / ((base_hex_coding) * (base_hex_coding)))

#ifndef DEBUG_ORDER
	#define LIMIT_HEX_ORDER(o)
#else
	#define LIMIT_HEX_ORDER(o) 						o = MAKE_HEX_ORDER(MAX_QUAD_ORDER, MAX_QUAD_ORDER, MAX_QUAD_ORDER);
#endif



inline Order1 get_hex_order_i(Order3 order, int i) {
	switch (i) {
		case 1: return GET_HEX_ORDER_1(order);
		case 2: return GET_HEX_ORDER_2(order);
		case 3: return GET_HEX_ORDER_3(order);
		default: assert(0); return -1;
	}
}

inline Order3 add_hex_orders(Order3 ord1, Order3 ord2) {
	int o1 = GET_HEX_ORDER_1(ord1) + GET_HEX_ORDER_1(ord2);
	if (o1 > MAX_QUAD_ORDER) o1 = MAX_QUAD_ORDER;
	int o2 = GET_HEX_ORDER_2(ord1) + GET_HEX_ORDER_2(ord2);
	if (o2 > MAX_QUAD_ORDER) o2 = MAX_QUAD_ORDER;
	int o3 = GET_HEX_ORDER_3(ord1) + GET_HEX_ORDER_3(ord2);
	if (o3 > MAX_QUAD_ORDER) o3 = MAX_QUAD_ORDER;
	return MAKE_HEX_ORDER(o1, o2, o3);
}

// Multiply orders in each direction by constant
inline Order3 mul_hex_orders(Order3 ord, int c) {
	int o1 = c * GET_HEX_ORDER_1(ord);
	if (o1 > MAX_QUAD_ORDER) o1 = MAX_QUAD_ORDER;
	int o2 = c * GET_HEX_ORDER_2(ord);
	if (o2 > MAX_QUAD_ORDER) o2 = MAX_QUAD_ORDER;
	int o3 = c * GET_HEX_ORDER_3(ord);
	if (o3 > MAX_QUAD_ORDER) o3 = MAX_QUAD_ORDER;
	return MAKE_HEX_ORDER(o1, o2, o3);
}

inline Order3 max_hex_order(Order3 ord1, Order3 ord2) {
	int o1 = GET_HEX_ORDER_1(ord1);
	if (o1 < GET_HEX_ORDER_1(ord2)) o1 = GET_HEX_ORDER_1(ord2);
	int o2 = GET_HEX_ORDER_2(ord1);
	if (o2 < GET_HEX_ORDER_2(ord2)) o2 = GET_HEX_ORDER_2(ord2);
	int o3 = GET_HEX_ORDER_3(ord1);
	if (o3 < GET_HEX_ORDER_3(ord2)) o3 = GET_HEX_ORDER_3(ord2);
	return MAKE_HEX_ORDER(o1, o2, o3);
}

inline int get_hex_face_order(int face, Order3 order) {
	Order1 o1 = GET_HEX_ORDER_1(order);
	Order1 o2 = GET_HEX_ORDER_2(order);
	Order1 o3 = GET_HEX_ORDER_3(order);

	if ((face == 0) || (face == 1))
		return MAKE_QUAD_ORDER(o2, o3);
	else if ((face == 2) || (face == 3))
		return MAKE_QUAD_ORDER(o1, o3);
	else if ((face == 4) || (face == 5))
		return MAKE_QUAD_ORDER(o1, o2);
	else
		EXIT(ERR_FACE_INDEX_OUT_OF_RANGE);
}

inline int get_hex_edge_order(int edge, Order3 order) {
	if ((edge == 0) || (edge == 2) || (edge == 10) || (edge == 8))
		return GET_HEX_ORDER_1(order);
	else if((edge == 1) || (edge == 3) || (edge == 11) || (edge == 9))
		return GET_HEX_ORDER_2(order);
	else if((edge == 4) || (edge == 5) || (edge == 6) || (edge == 7))
		return GET_HEX_ORDER_3(order);
	else
		EXIT(ERR_EDGE_INDEX_OUT_OF_RANGE);
}

inline Order3 turn_hex_face_order(Order3 ord) {
	Order1 o1 = GET_HEX_ORDER_1(ord);
	Order1 o2 = GET_HEX_ORDER_2(ord);
	Order1 o3 = GET_HEX_ORDER_3(ord);
	if (o1 <= 1) swapint(o2, o3);
	if (o2 <= 1) swapint(o1, o3);
	if (o3 <= 1) swapint(o1, o2);
	return MAKE_HEX_ORDER(o1, o2, o3);
}

// use izotropic order
inline Order3 make_hex_iso(Order3 order) {
	Order1 o1 = GET_HEX_ORDER_1(order);
	Order1 o2 = GET_HEX_ORDER_2(order);
	Order1 o3 = GET_HEX_ORDER_3(order);
	Order1 max = (o1 > o2 ? (o1 > o3 ? o1 : o3) : (o2 > o3 ? o2 : o3));
	return MAKE_HEX_ORDER(max, max, max);
}


// PRISM specific /////////////////////////////////////////////////////////////////////////////////

const int prism_base_coding = MAX_QUAD_ORDER_TRI + 1;

#define MAKE_PRISM_ORDER(h_order, v_order) ((v_order) * (prism_base_coding) + (h_order))
#define GET_PRISM_ORDER_1(order) ((order) % (prism_base_coding))
#define GET_PRISM_ORDER_2(order) ((order) / (prism_base_coding))

// TODO: prism specific macros/functions for faces, etc.

// TODO: limit prism order
#define LIMIT_PRISM_ORDER(o)

// Surface order //////////////////////////////////////////////////////////////////////////////////

enum EOrderType {
	OT_ELEM = 0,
	OT_FACE = 1,
	OT_EDGE = 2,
	OT_VERTEX = 3
};

//  Order: |--|----|--- .. ---|
//                 ++++++++++++ -- order
//            +++++             -- face/edge num
//          ++                  -- order type (EorderType)
//
#define MAKE_VOL_ORDER(order) ((order))
#define MAKE_FACE_ORDER(face, order) ((OT_FACE << 30) | (face << 26) | (order))
#define MAKE_EDGE_ORDER(edge, order) ((OT_EDGE << 30) | (edge << 26) | (order))
#define MAKE_VERTEX_ORDER() ((OT_VERTEX << 30))

#define GET_ORDER_TYPE(order) (EOrderType) ((order >> 30) & 0x03)
#define GET_EDGE_FROM_ORDER(order) ((order >> 26) & 0x0F)
#define GET_FACE_FROM_ORDER(order) ((order >> 26) & 0x0F)
#define GET_ORDER_FROM_ORDER(order) (order & 0x1FFFFFF)


#ifndef DEBUG_ORDER
	#define LIMIT_TRI_ORDER(o)
	#define LIMIT_QUAD_ORDER(o)
#else
	#define LIMIT_TRI_ORDER(o) 							o = MAX_QUAD_ORDER_TRI;
	#define LIMIT_QUAD_ORDER(o) 						o = MAKE_QUAD_ORDER(MAX_QUAD_ORDER, MAX_QUAD_ORDER);
#endif

#endif
