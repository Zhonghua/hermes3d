#ifndef _SHAPESET_H1_LOBATTO_TETRA_H_
#define _SHAPESET_H1_LOBATTO_TETRA_H_

#include "../shapeset.h"
#include "../refdomain.h"
#include "tetra.h"


/// H1 Lobatto shapeset for tetrahedra
///
/// @ingroup shapesets
class H1ShapesetLobattoTetra : public Shapeset {
public:
	H1ShapesetLobattoTetra();
	virtual ~H1ShapesetLobattoTetra();

	virtual int get_vertex_index(int vertex) const {
		CHECK_VERTEX(vertex);
		return vertex_indices[vertex];
	}

	virtual int *get_edge_indices(int edge, int ori, order1_t order) {
		CHECK_EDGE(edge); CHECK_ORDER(order);
		return edge_indices[edge][ori];
	}

	virtual int *get_face_indices(int face, int ori, order2_t order) {
		return face_indices[face][ori];
	}

	virtual int *get_bubble_indices(order3_t order) {
		return bubble_indices[order.get_idx()];
	}

	virtual int get_num_edge_fns(order1_t order) const {
		return edge_count[order];
	}

	virtual int get_num_face_fns(order2_t order) const {
		return face_count[order.get_idx()];
	}

	virtual int get_num_bubble_fns(order3_t order) const {
		return bubble_count[order.get_idx()];
	}

	virtual int get_face_orientations(int face) const { return RefTetra::get_face_orientations(face); }

	virtual int get_edge_orientations() const { return RefTetra::get_edge_orientations(); }

	virtual order3_t get_order(int index) const;

	virtual int get_shape_type(int index) const {
		return -1;
	}

protected:
	shape_fn_t **shape_table[VALUE_TYPES];

	/// Indices of vertex shape functions on reference element, indexing: []
	int   *vertex_indices;
	/// Indices of edge shape functions on reference element, indexing: [edge index][ori][]
	int ***edge_indices;
	/// Indices of face shape functions on reference element, indexing: [face index][ori][]
	int ***face_indices;
	/// Indices of bubble functions on reference element, indexing: [order][]
	int  **bubble_indices;

	/// Number of edge shape functions on reference element, indexing: [order]
	int *edge_count;
	/// Number of face shape functions on reference element, indexing: [order]
	int *face_count;
	/// Number of bubble functions on reference element, indexing: [order]
	int *bubble_count;

	virtual double get_val(int n, int index, double x, double y, double z, int component) {
		CHECK_COMPONENT(component);
		return shape_table[n][component][index](x, y, z);
	}
};

#undef CHECK_VERTEX
#undef CHECK_EDGE
#undef CHECK_FACE
#undef CHECK_FACE_MODE
#undef CHECK_FACE_ORI

#endif
