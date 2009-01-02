#ifndef _REFMAP_SHAPESET__H_
#define _REFMAP_SHAPESET__H_

/// Special shapeset used in RefMap for tetrahedra
///
/// @ingroup shapesets
class RefMapShapesetTetra : public Shapeset {
public:
	RefMapShapesetTetra();
	virtual ~RefMapShapesetTetra();

	// @return index of a vertex shape function for a vertex
	// @param [in] vertex - index of the vertex
	virtual int get_vertex_index(int vertex) const {
		return vertex_indices[vertex];
	}

	/// @return indices of edge shape functions
	/// @param [in] edge - edge number (local)
	/// @param [in] ori - orientation of the edge (0 or 1)
	/// @param [in] order - order on the edge
	virtual int *get_edge_indices(int edge, int ori, order1_t order) {
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	/// @return indices of face shape functions
	/// @param [in] face - face number (local)
	/// @param [in] ori - orinetation of the face
	/// @param [in] order - order on the face
	virtual int *get_face_indices(int face, int ori, order2_t order) {
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	/// @return indices of bubble functions
	/// @param order - order of the bubble function
	virtual int *get_bubble_indices(order3_t order) {
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	virtual int get_num_edge_fns(order1_t order) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual int get_num_face_fns(order2_t order) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual int get_num_bubble_fns(order3_t order) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual int get_face_orientations(int face) const { return 0; }

	virtual int get_edge_orientations() const { return 0; }

	virtual order3_t get_order(int index) const {
		CHECK_INDEX(index);
		return index_to_order[index];
	}

	virtual int get_shape_type(int index) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual double get_value(int n, int index, double x, double y, double z, int component) {
		CHECK_INDEX(index); CHECK_COMPONENT(component);
		return shape_table[n][component][index](x, y, z);
	}

protected:
	shape_fn_t **shape_table[VALUE_TYPES];

	/// Indices of vertex shape functions on reference element, indexing: [vertex shape fn index]
	int  *vertex_indices;

	virtual double get_val(int n, int index, double x, double y, double z, int component) {
		CHECK_INDEX(index); CHECK_COMPONENT(component);
		return shape_table[n][component][index](x, y, z);
	}
};


/// Special shapeset used in RefMap for hexahedra
///
/// @ingroup shapesets
class RefMapShapesetHex : public Shapeset {
public:
	RefMapShapesetHex();
	virtual ~RefMapShapesetHex();

	// @return index of a vertex shape function for a vertex
	// @param [in] vertex - index of the vertex
	virtual int get_vertex_index(int vertex) const {
		return vertex_indices[vertex];
	}

	/// @return indices of edge shape functions
	/// @param [in] edge - edge number (local)
	/// @param [in] ori - orientation of the edge (0 or 1)
	/// @param [in] order - order on the edge
	virtual int *get_edge_indices(int edge, int ori, order1_t order) {
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	/// @return indices of face shape functions
	/// @param [in] face - face number (local)
	/// @param [in] ori - orinetation of the face
	/// @param [in] order - order on the face
	virtual int *get_face_indices(int face, int ori, order2_t order) {
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	/// @return indices of bubble functions
	/// @param order - order of the bubble function
	virtual int *get_bubble_indices(order3_t order) {
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	virtual int get_num_edge_fns(order1_t order) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual int get_num_face_fns(order2_t order) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual int get_num_bubble_fns(order3_t order) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual int get_face_orientations(int face) const { return 0; }

	virtual int get_edge_orientations() const { return 0; }

	virtual order3_t get_order(int index) const {
		CHECK_INDEX(index);
		return order;
	}

	virtual int get_shape_type(int index) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual double get_value(int n, int index, double x, double y, double z, int component) {
		CHECK_INDEX(index); CHECK_COMPONENT(component);
		return shape_table[n][component][index](x, y, z);
	}

protected:
	shape_fn_t **shape_table[VALUE_TYPES];
	order3_t order;

	/// Indices of vertex shape functions on reference element, indexing: [vertex shape fn index]
	int  *vertex_indices;

	virtual double get_val(int n, int index, double x, double y, double z, int component) {
		CHECK_INDEX(index); CHECK_COMPONENT(component);
		return shape_table[n][component][index](x, y, z);
	}
};

#undef CHECK_VERTEX

#endif
