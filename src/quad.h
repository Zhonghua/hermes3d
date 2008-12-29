#ifndef _QUAD_H_
#define _QUAD_H_

#include "config.h"
#include "common.h"
#include "order.h"
#include <common/error.h>
#include <common/array.h>

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
};


//
// 3D quadratures
//


/// Numerical quadratures in 3D
///
/// @ingroup quadratures
class Quad3D {
public:
	virtual QuadPt3D *get_points(order3_t order) { return tables[order.get_idx()]; }
	virtual int get_num_points(order3_t order) { return np[order.get_idx()]; }

	virtual QuadPt3D *get_edge_points(int edge, int order) { return edge_tables[edge][order]; }
	int get_edge_num_points(int order) const { return np_edge[order]; }

	virtual QuadPt3D *get_face_points(int face, order2_t order) { return face_tables[face][order.get_idx()]; }
	int get_face_num_points(int face, order2_t order) const { return np_face[order.get_idx()]; }

	virtual QuadPt3D *get_vertex_points() { return vertex_table; }
	int get_vertex_num_points() const { return np_vertex; }

	order3_t get_max_order() const { return max_order; }
	int get_edge_max_order(int edge) const { return max_edge_order; }
	order2_t get_face_max_order(int face) const { return max_face_order; }

	EMode3D get_mode() const { return mode; }

protected:
	/// mode of quadratures (MODE_TETRAHEDRON, MODE_HEXAHEDRON, MODE_PRISM)
	EMode3D mode;
	/// maximal order for integration (interpretation depeneds on the mode)
	order3_t max_order;
	int max_edge_order;
	order2_t max_face_order;

	Array<QuadPt3D *> tables;
	QuadPt3D ***edge_tables;
	QuadPt3D ***face_tables;
	QuadPt3D *vertex_table;
	Array<int> np;
	int *np_edge;
	int *np_face;
	int np_vertex;
};


// interface for getting quadratures - library wide ////////////////////////////////////////////////

Quad1D *get_quadrature(EMode1D mode);
Quad2D *get_quadrature(EMode2D mode);
Quad3D *get_quadrature(EMode3D mode);


// Helpers for calculating with orders ////////////////////////////////////////////////////////////


// QUAD specific //////////////////////////////////////////////////////////////////////////////////

//const int base_quad_coding = MAX_QUAD_ORDER + 1;

//#define MAKE_QUAD_ORDER(h_order, v_order) ((v_order) * (base_quad_coding) + (h_order))
//#define GET_QUAD_ORDER_1(order) ((order) % (base_quad_coding))
//#define GET_QUAD_ORDER_2(order) ((order) / (base_quad_coding))

//inline Order2 add_quad_orders(Order2 ord1, Order2 ord2) {
//	int o1 = GET_QUAD_ORDER_1(ord1) + GET_QUAD_ORDER_1(ord2);
//	if (o1 > MAX_QUAD_ORDER) o1 = MAX_QUAD_ORDER;
//	int o2 = GET_QUAD_ORDER_2(ord1) + GET_QUAD_ORDER_2(ord2);
//	if (o2 > MAX_QUAD_ORDER) o2 = MAX_QUAD_ORDER;
//	return MAKE_QUAD_ORDER(o1, o2);
//}

//inline void make_quad_debug(Order2 &order) {
//	Order1 o1 = GET_QUAD_ORDER_1(order);
//	Order1 o2 = GET_QUAD_ORDER_2(order);
//	Order1 max = o1 > o2 ? o1 : o2;
//	order = MAKE_QUAD_ORDER(max, max);
//}

// TETRA specific /////////////////////////////////////////////////////////////////////////////////

//#ifndef DEBUG_ORDER
//	#define LIMIT_TETRA_ORDER(o) 						if ((o) > MAX_QUAD_ORDER_TETRA) o = MAX_QUAD_ORDER_TETRA;
//#else
//	#define LIMIT_TETRA_ORDER(o) 						o = MAX_QUAD_ORDER_TETRA;
//#endif


// HEX specific ///////////////////////////////////////////////////////////////////////////////////

//const int base_hex_coding = MAX_QUAD_ORDER + 1;

//#define MAKE_HEX_ORDER(h_order, f_order, v_order) (((v_order) * (base_hex_coding) + (f_order)) * (base_hex_coding) + (h_order))
//#define GET_HEX_ORDER_1(order) ((order) % (base_hex_coding))
//#define GET_HEX_ORDER_2(order) (((order) / (base_hex_coding)) % (base_hex_coding))
//#define GET_HEX_ORDER_3(order) ((order) / ((base_hex_coding) * (base_hex_coding)))

//#ifndef DEBUG_ORDER
//	#define LIMIT_HEX_ORDER(o)
//#else
//	#define LIMIT_HEX_ORDER(o) 						o = MAKE_HEX_ORDER(MAX_QUAD_ORDER, MAX_QUAD_ORDER, MAX_QUAD_ORDER);
//#endif


/*
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
*/

// PRISM specific /////////////////////////////////////////////////////////////////////////////////

//const int prism_base_coding = MAX_QUAD_ORDER_TRI + 1;
//
//#define MAKE_PRISM_ORDER(h_order, v_order) ((v_order) * (prism_base_coding) + (h_order))
//#define GET_PRISM_ORDER_1(order) ((order) % (prism_base_coding))
//#define GET_PRISM_ORDER_2(order) ((order) / (prism_base_coding))

// TODO: prism specific macros/functions for faces, etc.

// TODO: limit prism order
#define LIMIT_PRISM_ORDER(o)

// Surface order //////////////////////////////////////////////////////////////////////////////////

//// Qorder ////
//
// quadrature order
//

#define QOT_ELEMENT						0
#define QOT_FACE						1
#define QOT_EDGE						2
#define QOT_VERTEX						3

struct qorder_t {
	unsigned type:3;				// QOT_XXX
	union {
		unsigned edge: 4;			// the number of the local edge (if type == QOT_EDGE)
		unsigned face: 4;			// the number of the local face (if type == QOT_FACE)
	};
	unsigned order:25;				// order (TODO: )

	qorder_t(unsigned type) {
		this->type = type;
	}

	qorder_t(unsigned type, order3_t order) {
		this->type = type;
//		this->order = order;
	}

	qorder_t(unsigned type, unsigned ef, order2_t order) {
		this->type = type;
		this->edge = ef;
//		this->order = order;
	}

	operator int() { return (((type << 4) | edge) << 25) | order; }
};

#define ELEM_QORDER(o)				qorder_t(QOT_ELEMENT, o)
#define FACE_QORDER(f, o)			qorder_t(QOT_FACE, f, o)
#define EDGE_QORDER(e, o)			qorder_t(QOT_EDGE, e, o)
#define VTX_QORDER()				qorder_t(QOT_VERTEX)

//

#ifndef DEBUG_ORDER
	#define LIMIT_TRI_ORDER(o)
	#define LIMIT_QUAD_ORDER(o)
#else
	#define LIMIT_TRI_ORDER(o) 							o = MAX_QUAD_ORDER_TRI;
	#define LIMIT_QUAD_ORDER(o) 						o = MAKE_QUAD_ORDER(MAX_QUAD_ORDER, MAX_QUAD_ORDER);
#endif

#endif
