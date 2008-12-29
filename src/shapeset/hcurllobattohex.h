#ifndef _HCURL_SHAPESET_LOBATTO_HEX_H_
#define _HCURL_SHAPESET_LOBATTO_HEX_H_

#include "../shapeset.h"
#include "../mesh.h"
#include "hex.h"

#define FN_TYPE_EDGE 0
#define FN_TYPE_FACE 1
#define FN_TYPE_BUBBLE 2

///
/// @brief H1 shapeset for hexahedra
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
class HCurlShapesetLobattoHex : public Shapeset {
public:
	HCurlShapesetLobattoHex();
	virtual ~HCurlShapesetLobattoHex();

	virtual int get_vertex_index(int vertex) const {
		assert(false); //no vertex functions in hcurl
	}

	virtual int *get_edge_indices(int edge, int ori, int order) {
		CHECK_EDGE(edge); CHECK_EDGE_ORDER(order);
//		if (order == -1) order = max_edge_order;
		if (!edge_indices[edge][ori].exists(order)) compute_edge_indices(edge, ori, order);
		return edge_indices[edge][ori][order];
	}

	virtual int *get_face_indices(int face, int ori, order2_t order) {
//		CHECK_FACE(face); CHECK_FACE_ORDER(order);
//		if (order == -1) order = max_face_order;
		if (!face_indices[face][ori].exists(order.get_idx())) compute_face_indices(face, ori, order);
		return face_indices[face][ori][order.get_idx()];
	}

	virtual int *get_bubble_indices(order3_t order) {
//		CHECK_ORDER(order);
//		if (order == -1) order = max_order;
		if (!bubble_indices.exists(order.get_idx())) compute_bubble_indices(order);
		return bubble_indices[order.get_idx()];
	}

	virtual int get_num_edge_fns(int order) const {
		CHECK_EDGE_ORDER(order);
		return (order + 1);
	}

	virtual int get_num_face_fns(order2_t order) const {
//		CHECK_FACE_ORDER(order);
		int order1 = order.x;
		int order2 = order.y;

		return (order1 + 1) * order2 + order1 * (order2 + 1);
	}

	virtual int get_num_bubble_fns(order3_t order) const {
//		CHECK_ORDER(order);
		int order1 = order.x;
		int order2 = order.y;
		int order3 = order.z;

		return (order1 + 1) * order2 * order3 + order1 * (order2 + 1) * order3 + order1 * order2 * (order3 + 1);
	}

	virtual int get_face_orientations(int face) const { return RefHex::get_face_orientations(face); }

	virtual int get_edge_orientations() const { return RefHex::get_edge_orientations(); }

	virtual order3_t get_order(int index) const;

	virtual int get_shape_type(int index) const {
		return -1;
	}

	// returns 0 if vectors in direction with smaller number of the face, 1 otherwise
	virtual int get_facefn_variant(int index) const {
		int idx = GET_IDX_FROM_INDEX(index);
		int ind[3], fn_type, unit_index, which_legendre, v_direction;
		check_fn_index(index, ind, fn_type, unit_index, which_legendre, v_direction);
		assert(fn_type == SHAPE_FN_FACE);
		switch(unit_index){
			case 0 :
			case 1 : {
				if(v_direction == 1)
					return 0;
				else if(v_direction == 2)
					return 1;
				else
					assert(0);
			}
			case 2 :
			case 3 : {
				if(v_direction == 0)
					return 0;
				else if(v_direction == 2)
					return 1;
				else
					assert(0);
			}
			case 4 :
			case 5 : {
				if(v_direction == 0)
					return 0;
				else if(v_direction == 1)
					return 1;
				else
					assert(0);
			}

		}
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
	void compute_face_indices(int face, int ori, order2_t order);
	void compute_bubble_indices(order3_t order);

/*	virtual double get_val(int n, int index, double x, double y, double z, int component) {
		// use on-the-fly function
		if (shape_table_deleg[n] == NULL) EXIT(ERR_FAILURE, "Missing a delegate function for calculating shape functions");
		return shape_table_deleg[n](index, x, y, z, component);
	}
	/// --- put CED specific stuff here ---
	virtual CEDComb *calc_constrained_edge_combination(int ori, int order, Part part);
	virtual CEDComb *calc_constrained_edge_face_combination(int ori, int order, Part part);
	virtual CEDComb *calc_constrained_face_combination(int ori, int order, Part part);
*/

	public: //TODO protected
//protected:
	///  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	///TODO
	/// the reason, why all following functions are static is this:
	/// pointers to calc_fn_value etc. are used to define ..._deleg tables
	/// compiler would not allowe to do it, if it is not static member function (I dont know how)
	/// In H1 those functions are not member at all, but I need it here - using other member data.
	///  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	static void check_edge_fn_index(int edge_fn_index, int indices[3], int &edge_index, int &which_legendre, int &t_direction, int &v_direction);
	static void check_face_fn_index(int face_fn_index, int indices[3],  int &face_index,
									 int &which_legendre, int &t_direction_1, int &t_direction_2, int &v_direction);
	static void check_bubble_fn_index(int bubble_fn_index, int indices[3], int &which_legendre, int &v_direction);
	static void check_fn_index(int index, int indices[3], int &fn_type, int &unit_index,
							  int &which_legendre, int &v_direction);
	static void dump_fn_index(int index);

	///value(which_der=0) and partial derivatives(which_der=1,2,3) calculated in one function
	static double calc_any_value(int index, double x, double y, double z, int component, int which_der);

	static double calc_fn_value(int index, double x, double y, double z, int component);
	static double calc_dx_value(int index, double x, double y, double z, int component);
	static double calc_dy_value(int index, double x, double y, double z, int component);
	static double calc_dz_value(int index, double x, double y, double z, int component);

	static int max_fns_per_edge;
	static int max_fns_per_face;
	static int max_bubble_fns;

	static int face_offset;
	static int bubble_offset;

	virtual double get_val(int n, int index, double x, double y, double z, int component) {
		//		int ori = GET_ORI_FROM_INDEX(index);
		//		int idx = GET_IDX_FROM_INDEX(index);
		//		CHECK_INDEX(idx);
		// use on-the-fly function
		if (shape_table_deleg[n] == NULL) EXIT(ERR_FAILURE, "Missing a delegate function for calculating shape functions");
		return shape_table_deleg[n](index, x, y, z, component);
	}


// implementation has to be redefined from standart
// for face functions we use only part of facefunctions defined for given order
	virtual double get_constrained_value(int n, int index, double x, double y, double z, int component);
	double get_constrained_face_value(int n, int index, double x, double y, double z, int component);

	/// --- put CED specific stuff here ---
	virtual CEDComb *calc_constrained_edge_combination(int ori, int order, Part part);
	// facenf_type - 0 for vectors in first direction, 1 for second
	virtual CEDComb *calc_constrained_edge_face_combination(int ori, int order, Part part, int facefn_variant);
	virtual CEDComb *calc_constrained_face_combination(int ori, int order, Part part, int facefn_variant);
};

#undef CHECK_VERTEX
#undef CHECK_EDGE
#undef CHECK_FACE
#undef CHECK_FACE_MODE
#undef CHECK_FACE_ORI

#endif

