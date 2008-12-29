#ifndef _REFMAP_H_
#define _REFMAP_H_

#include "common.h"
#include "transform.h"
#include "quad.h"
#include "mesh.h"
#include "precalc.h"

/// @defgroup refmap Reference mapping
///
/// TODO: description


/// Reference mapping (for evaluating integrals on a physical domain)
///
/// @ingroup refmap
class RefMap : public Transformable {
public:
	/// NOTE: Use this constructor only if constructing an array of RefMaps, then call set_mesh immediately!
	RefMap();
	/// NOTE: Always use this constructor whenever it is possible
	RefMap(Mesh *mesh);
	virtual ~RefMap();

	/// Sets the mesh where the reference map is working
	/// @param mesh [in] Pointer to the mesh.
	void set_mesh(Mesh *mesh) { this->mesh = mesh; }

	/// Sets the quadrature points in which the reference map will be evaluated.
	/// @param quad [in] The quadrature points.
	void set_quad(Quad3D *quad);

	/// Initializes the reference map for the specified element.
	/// Must be called prior to using all other functions in the class.
	virtual void set_active_element(Element *e);

	/// \return True if the jacobian of the reference map is constant (which
	/// is the case for non-curvilinear tetra elements), false otherwise.
	bool is_jacobian_const() const { return is_const; }

//	/// \return True if the jacobian of the reference map is diagonal (even
//	/// though probably non-constant). This is the case for axis-aligned
//	/// non-curvilinear hex elements (???).
//	bool is_jacobian_diag() const { return is_diag; }

	/// \return The increase in the integration order due to the reference map.
	order3_t get_ref_order() const { return ref_order; }

	/// \return The increase in the integration order due to the inverse reference map.
	order3_t get_inv_ref_order() const { return inv_ref_order; }

	/// If the jacobian of the reference map is constant, this is the fast
	/// way to obtain it.
	/// @return The constant jacobian of the reference map.
	double get_const_jacobian() const { return const_jacobian; }

	/// If the reference map is constant, this is the fast way to obtain
	/// its jacobi matrix.
	/// @return The matrix of the reference map.
	double3x3 *get_const_ref_map() { return &const_ref_map; }

	/// If the reference map is constant, this is the fast way to obtain
	/// its inverse matrix.
	/// @return The transposed inverse matrix of the reference map.
	double3x3 *get_const_inv_ref_map() { return &const_inv_ref_map; }

	/// \return The jacobian of the reference map precalculated at the integration
	/// points of the specified order. Intended for non-constant jacobian elements.
	/// \param order [in] Integration order
	double *get_jacobian(qorder_t order) {
		if (!cur_node->jacobian.exists(order)) calc_inv_ref_map(order);
		return cur_node->jacobian[order];
	}

	/// \return The jacobi matrices of the reference map precalculated at the
	/// integration points of the specified order. Intended for non-constant
	/// jacobian elements.
	/// \param order [in] Integration order
	double3x3 *get_ref_map(qorder_t order) {
		if (!cur_node->inv_ref_map.exists(order)) calc_inv_ref_map(order);
		return cur_node->ref_map[order];
	}

	/// \return The transposed inverse matrices of the reference map precalculated at the
	/// integration points of the specified order. Intended for non-constant
	/// jacobian elements.
	/// \param order [in] Integration order
	double3x3 *get_inv_ref_map(qorder_t order) {
		assert(cur_node != NULL);
		if (!cur_node->inv_ref_map.exists(order)) calc_inv_ref_map(order);
		return cur_node->inv_ref_map[order];
	}

	/// \return The x-coordinates of the integration points transformed to the
	/// physical domain of the element. Intended for integrals containing spatial
	/// variables.
	/// \param order [in] Integration order
	double *get_phys_x(order3_t order) {
		if (!cur_node->phys_x.exists(order.get_idx())) calc_phys_x(order);
		return cur_node->phys_x.get(order.get_idx());
	}

	/// \return The y-coordinates of the integration points transformed to the
	/// physical domain of the element. Intended for integrals containing spatial
	/// variables.
	/// \param order [in] Integration order
	double *get_phys_y(order3_t order) {
		if (!cur_node->phys_y.exists(order.get_idx())) calc_phys_y(order);
		return cur_node->phys_y.get(order.get_idx());
	}

	/// \return The z-coordinates of the integration points transformed to the
	/// physical domain of the element. Intended for integrals containing spatial
	/// variables.
	/// \param order [in] Integration order
	double *get_phys_z(order3_t order) {
		if (!cur_node->phys_z.exists(order.get_idx())) calc_phys_z(order);
		return cur_node->phys_z.get(order.get_idx());
	}

	// Edges //////////////////////////////////////////////////////////////////////////////////////

	/// @return The x-coordinates of the integration points transformed to the
	/// physical domain of the element. Intended for integrals containing spatial
	/// variables.
	/// @param edge [in] Number of edge of the element
	/// @param order [in] Integration order
	double *get_edge_phys_x(int edge, order1_t order) {
		if (!cur_node->edge_phys_x[edge].exists(order)) calc_edge_phys_x(edge, order);
		return cur_node->edge_phys_x[edge].get(order);
	}

	/// @return The y-coordinates of the integration points transformed to the
	/// physical domain of the element. Intended for integrals containing spatial
	/// variables.
	/// @param edge [in] Number of edge of the element
	/// @param order [in] Integration order
	double *get_edge_phys_y(int edge, order1_t order) {
		if (!cur_node->edge_phys_y[edge].exists(order)) calc_edge_phys_y(edge, order);
		return cur_node->edge_phys_y[edge].get(order);
	}

	/// @return The z-coordinates of the integration points transformed to the
	/// physical domain of the element. Intended for integrals containing spatial
	/// variables.
	/// @param edge [in] Number of edge of the element
	/// @param order [in] Integration order
	double *get_edge_phys_z(int edge, order1_t order) {
		if (!cur_node->edge_phys_z[edge].exists(order)) calc_edge_phys_z(edge, order);
		return cur_node->edge_phys_z[edge].get(order);
	}

	// Faces //////////////////////////////////////////////////////////////////////////////////////

	inline bool is_face_const_jacobian(int face) const {
		return (cur_node->face_mode[face] == MODE_TRIANGLE);
	}

	double get_face_const_jacobian(int face) {
		if (cur_node->face_const_jacobian[face] == 0) calc_face_const_jacobian(face);
		return cur_node->face_const_jacobian[face];
	}

	double *get_face_jacobian(int face, order2_t order) {
		if (!cur_node->face_jacobian[face].exists(order.get_idx())) calc_face_jacobian(face, order);
		return cur_node->face_jacobian[face].get(order.get_idx());
	}

	double3x3 *get_face_const_ref_map(int face) {
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	double3x3 *get_face_ref_map(int face, order2_t order) {
		if (!cur_node->face_ref_map[face].exists(order.get_idx())) calc_face_inv_ref_map(face, order);
		return cur_node->face_ref_map[face][order.get_idx()];
	}

	double3x3 *get_face_const_inv_ref_map(int face){
		EXIT(ERR_NOT_IMPLEMENTED);
		return NULL;
	}

	double3x3 *get_face_inv_ref_map(int face, order2_t order) {
		if (!cur_node->face_inv_ref_map[face].exists(order.get_idx())) calc_face_inv_ref_map(face, order);
		return cur_node->face_inv_ref_map[face][order.get_idx()];
	}

	// outer unit normal to the face (transformed)
 	Point3D *get_face_const_normal(int face) {
		EXIT(ERR_NOT_IMPLEMENTED);
 		//if (face_const_normal[face] == 0) calc_face_const_normal(face);
 		return &(cur_node->face_const_normal[face]);
 	}

	//outer unit normal to the face (transformed)
 	Point3D *get_face_normal(int face, order2_t order) {
 		if (!cur_node->face_normal[face].exists(order.get_idx())) calc_face_normal(face, order);
 		return cur_node->face_normal[face].get(order.get_idx());
 	}

	double *get_face_phys_x(int face, order2_t order) {
		if (!cur_node->face_phys_x[face].exists(order.get_idx())) calc_face_phys_x(face, order);
		return cur_node->face_phys_x[face].get(order.get_idx());
	}

	double *get_face_phys_y(int face, order2_t order) {
		if (!cur_node->face_phys_y[face].exists(order.get_idx())) calc_face_phys_y(face, order);
		return cur_node->face_phys_y[face].get(order.get_idx());
	}

	double *get_face_phys_z(int face, order2_t order) {
		if (!cur_node->face_phys_z[face].exists(order.get_idx())) calc_face_phys_z(face, order);
		return cur_node->face_phys_z[face].get(order.get_idx());
	}

	// Vertices /////////

	double *get_vertex_phys_x() {
		if (cur_node->vertex_phys_x == NULL) calc_vertex_phys();
		return cur_node->vertex_phys_x;
	}

	double *get_vertex_phys_y() {
		if (cur_node->vertex_phys_y == NULL) calc_vertex_phys();
		return cur_node->vertex_phys_y;
	}

	double *get_vertex_phys_z() {
		if (cur_node->vertex_phys_z == NULL) calc_vertex_phys();
		return cur_node->vertex_phys_z;
	}

	// Transform ////////

	/// See Transformable::push_transform()
	virtual void push_transform(int son);

	/// See Transformable::pop_transform()
	virtual void pop_transform();

	void force_transform(uint64 sub_idx, Trf *ctm);

	void free();

protected:
	Mesh *mesh;
	Quad3D *quad;

	PrecalcShapeset *pss;

	bool      is_const;
	double    const_jacobian;
	double3x3 const_inv_ref_map;
	double3x3 const_ref_map;

	order3_t ref_order;
	order3_t inv_ref_order;

	struct Node {
		ArrayPtr<double>    jacobian;

		ArrayPtr<double3x3> inv_ref_map;
		ArrayPtr<double3x3> ref_map;

		Array<double *>     phys_x;
		Array<double *>     phys_y;
		Array<double *>     phys_z;

		Array<double *>     *edge_phys_x;
		Array<double *>     *edge_phys_y;
		Array<double *>     *edge_phys_z;

		Array<double *>     *face_phys_x;
		Array<double *>     *face_phys_y;
		Array<double *>     *face_phys_z;
		EMode2D             *face_mode;
		ArrayPtr<double>    *face_jacobian;
		double              *face_const_jacobian;
		ArrayPtr<double3x3> *face_ref_map;
		double3x3           *face_const_ref_map;
		ArrayPtr<double3x3> *face_inv_ref_map;
		double3x3           *face_const_inv_ref_map;
		ArrayPtr<Point3D>   *face_normal;
		Point3D             *face_const_normal;

		double				*vertex_phys_x;
		double				*vertex_phys_y;
		double				*vertex_phys_z;

		// for internal use
		int nedges, nfaces;
	};

	void calc_inv_ref_map(qorder_t order);
	void calc_const_inv_ref_map();
	int  calc_inv_ref_order();

	void calc_phys_x(order3_t order);
	void calc_phys_y(order3_t order);
	void calc_phys_z(order3_t order);
//	void calc_tangent(int edge);

	void calc_edge_phys_x(int edge, order1_t order);
	void calc_edge_phys_y(int edge, order1_t order);
	void calc_edge_phys_z(int edge, order1_t order);

	void calc_face_const_jacobian(int face);
	void calc_face_jacobian(int face, order2_t order);
	void calc_face_normal(int face, order2_t order);
	void calc_face_inv_ref_map(int face, order2_t order);
	void calc_face_phys_x(int face, order2_t order);
	void calc_face_phys_y(int face, order2_t order);
	void calc_face_phys_z(int face, order2_t order);

	void calc_vertex_phys();

	void *nodes;
	Node *cur_node;
	Node *overflow;

	void init_node(Node **pp);
	void free_node(Node *node);
	Node **handle_overflow();
	void update_cur_node();

	int num_coefs;					// # of coeffs in 'indices' array
	int indices[70];				// FIXME: magic number

	Vertex *coefs;
	Vertex vertex[8];				// max number of vertices (quad has 8 vertices, other elements have less)

	friend class ExactSolution;
};

#endif
