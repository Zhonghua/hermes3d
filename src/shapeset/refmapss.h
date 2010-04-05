// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#ifndef _REFMAP_SHAPESET__H_
#define _REFMAP_SHAPESET__H_

#include "common.h"
#include "../shapeset.h"
#include "../mesh.h"
#include "../refdomain.h"
#include "hex.h"

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
		return order3_t(1);
	}

	virtual order3_t get_dcmp(int index) const { return order3_t(-1); }

	virtual int get_shape_type(int index) const {
		EXIT(ERR_NOT_IMPLEMENTED);
		return 0;
	}

	virtual void get_values(int n, int index, int np, QuadPt3D *pt, int component, double *vals) {
		CHECK_COMPONENT(component);
		for (int k = 0; k < np; k++)
			vals[k] = shape_table[n][component][index](pt[k].x, pt[k].y, pt[k].z);
	}

	virtual double get_value(int n, int index, double x, double y, double z, int component) {
		CHECK_COMPONENT(component);
		return shape_table[n][component][index](x, y, z);
	}

protected:
	shape_fn_t **shape_table[VALUE_TYPES];

	/// Indices of vertex shape functions on reference element, indexing: [vertex shape fn index]
	int  *vertex_indices;
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
		CHECK_EDGE(edge);
		if (!edge_indices[edge][ori].exists(order)) compute_edge_indices(edge, ori, order);
		return edge_indices[edge][ori][order];
	}

	/// @return indices of face shape functions
	/// @param [in] face - face number (local)
	/// @param [in] ori - orinetation of the face
	/// @param [in] order - order on the face
	virtual int *get_face_indices(int face, int ori, order2_t order) {
		CHECK_FACE(face);
		if (!face_indices[face][ori].exists(order.get_idx())) compute_face_indices(face, ori, order);
		return face_indices[face][ori][order.get_idx()];
	}

	/// @return indices of bubble functions
	/// @param order - order of the bubble function
  	virtual int *get_bubble_indices(order3_t order) {
		if (!bubble_indices.exists(order.get_idx())) compute_bubble_indices(order);
		return bubble_indices[order.get_idx()];
	}

	virtual int get_num_edge_fns(order1_t order) const { return (order > 1) ? order - 1 : 0; }

	virtual int get_num_face_fns(order2_t order) const { return (order.x > 1 && order.y > 1) ? (order.x - 1) * (order.y - 1) : 0; }

	virtual int get_num_bubble_fns(order3_t order) const {
		return (order.x > 1 && order.y > 1 && order.z > 1) ? (order.x - 1) * (order.y - 1) * (order.z - 1) : 0;
	}

	virtual int get_face_orientations(int face) const { return RefHex::get_face_orientations(face); }

	virtual int get_edge_orientations() const { return RefHex::get_edge_orientations(); }

	virtual order3_t get_order(int index) const;

	virtual order3_t get_dcmp(int index) const;

	virtual int get_shape_type(int index) const;

	virtual void get_values(int n, int index, int np, QuadPt3D *pt, int component, double *vals) {
		if (index >= 0) shape_table_deleg[n](index, np, pt, component, vals);
		else get_constrained_values(n, index, np, pt, component, vals);
	}

	virtual double get_value(int n, int index, double x, double y, double z, int component) {
		QuadPt3D one(x, y, z, 1.0);
		double val = 0.0;
		if (index >= 0) shape_table_deleg[n](index, 1, &one, component, &val);
		else val = get_constrained_value(n, index, x, y, z, component);
		return val;
	}

protected:
	// some constants
	static const int NUM_EDGE_ORIS = 2;
	static const int NUM_FACE_ORIS = 8;

	shape_fn_deleg_t shape_table_deleg[VALUE_TYPES];

	/// Indices of vertex shape functions on reference element, indexing: [vertex shape fn index]
	int *vertex_indices;
	Array<int *> edge_indices[Hex::NUM_EDGES][NUM_EDGE_ORIS];
	Array<int *> face_indices[Hex::NUM_FACES][NUM_FACE_ORIS];
	Array<int *> bubble_indices;

	void compute_edge_indices(int edge, int ori, order1_t order);
	void compute_face_indices(int face, int ori, order2_t order);
	void compute_bubble_indices(order3_t order);
};

#undef CHECK_COMPONENT

#endif
