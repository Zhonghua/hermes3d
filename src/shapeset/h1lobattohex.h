#ifndef _H1_SHAPESET_LOBATTO_HEX_H_
#define _H1_SHAPESET_LOBATTO_HEX_H_

#include "../shapeset.h"
#include "../mesh.h"
#include "../refdomain.h"
#include "hex.h"

/// H1 shapeset for hexahedra
///
/// NOTE: This shapeset is very large. We use product geometry for calculating values of shape functions. Each function index
/// has to represent a shape function and its orientation. The encoeing is as following:
///  iii...iiii|ooo
///
///  ooo - last 3 bits represent the orientation (there is 8 possible orientations)
///  iii - is the index of shape function
///
/// Use GET_ORI_FROM_INDEX and GET_IDX_FROM_INDEX macroes for working with shape function indices.
///
/// @ingroup shapesets
class H1ShapesetLobattoHex : public Shapeset {
public:
	H1ShapesetLobattoHex();
	virtual ~H1ShapesetLobattoHex();

	virtual int get_vertex_index(int vertex) const {
		CHECK_VERTEX(vertex);
		return vertex_indices[vertex];
	}

	virtual int *get_edge_indices(int edge, int ori, int order) {
		CHECK_EDGE(edge); CHECK_EDGE_ORDER(order);
		if (order == 0) order = max_edge_order;
		if (!edge_indices[edge][ori].exists(order)) compute_edge_indices(edge, ori, order);
		return edge_indices[edge][ori][order];
	}

	virtual int *get_face_indices(int face, int ori, int order) {
 		CHECK_FACE(face); CHECK_FACE_ORDER(order);
		if (order == 0) order = max_face_order;
		if (!face_indices[face][ori].exists(order)) compute_face_indices(face, ori, order);
		return face_indices[face][ori][order];
	}

  	virtual int *get_bubble_indices(int order) {
 		CHECK_ORDER(order);
		 if (order == 0) order = max_order;
		if (!bubble_indices.exists(order)) compute_bubble_indices(order);
		return bubble_indices[order];
	}

	virtual int get_num_edge_fns(int order) const {
		CHECK_EDGE_ORDER(order);
		if (order > 1) return (order - 1);
		else return 0;
	}

	virtual int get_num_face_fns(int order) const {
		CHECK_FACE_ORDER(order);
		int order1 = GET_QUAD_ORDER_1(order);
		int order2 = GET_QUAD_ORDER_2(order);

		if (order1 > 1 && order2 > 1) return (order1 - 1) * (order2 - 1);
		else return 0;
	}

	virtual int get_num_bubble_fns(int order) const {
		CHECK_ORDER(order);
		int order1 = GET_HEX_ORDER_1(order);
		int order2 = GET_HEX_ORDER_2(order);
		int order3 = GET_HEX_ORDER_3(order);

		if (order1 > 1 && order2 > 1 && order3 > 1) return (order1 - 1) * (order2 - 1) * (order3 - 1);
		else return 0;
	}

	virtual int get_face_orientations(int face) const { return RefHex::get_face_orientations(face); }

	virtual int get_edge_orientations() const { return RefHex::get_edge_orientations(); }

	virtual int get_order(int index) const;

	virtual int get_shape_type(int index) const {
		return -1;
	}

protected:
	// some constants
	static const int NUM_EDGE_ORIS = 2;
	static const int NUM_FACE_ORIS = 8;

	shape_fn_deleg_t shape_table_deleg[VALUE_TYPES];

	// for validation (set in the constructor)
	int max_edge_order;
	int max_face_order;

	/// Indices of vertex shape functions on reference element, indexing: [vertex shape fn index]
	int *vertex_indices;
	Array<int *> edge_indices[Hex::NUM_EDGES][NUM_EDGE_ORIS];
	Array<int *> face_indices[Hex::NUM_FACES][NUM_FACE_ORIS];
	Array<int *> bubble_indices;

	void compute_edge_indices(int edge, int ori, int order);
	void compute_face_indices(int face, int ori, int order);
	void compute_bubble_indices(int order);

	virtual double get_val(int n, int index, double x, double y, double z, int component) {
//		int ori = GET_ORI_FROM_INDEX(index);
//		int idx = GET_IDX_FROM_INDEX(index);

//		CHECK_INDEX(idx);
		// use on-the-fly function
		if (shape_table_deleg[n] == NULL) EXIT(ERR_FAILURE, "Missing a delegate function for calculating shape functions");
		return shape_table_deleg[n](index, x, y, z, component);
	}


	/// --- put CED specific stuff here ---
	virtual CEDComb *calc_constrained_edge_combination(int ori, int order, Part part);
	virtual CEDComb *calc_constrained_edge_face_combination(int ori, int order, Part part, int dir);
	virtual CEDComb *calc_constrained_face_combination(int ori, int order, Part part);
};

#undef CHECK_VERTEX
#undef CHECK_EDGE
#undef CHECK_FACE
#undef CHECK_FACE_MODE
#undef CHECK_FACE_ORI

#endif

