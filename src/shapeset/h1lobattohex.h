// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
// Copyright (c) 2007 - 2009 Pavel Kus <pavel.kus@gmail.com>
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

#ifndef _H1_SHAPESET_LOBATTO_HEX_H_
#define _H1_SHAPESET_LOBATTO_HEX_H_

#include "../shapeset.h"
#include "../mesh.h"
#include "../refdomain.h"
#include "hex.h"

/// H1 shapeset for hexahedra
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

	virtual int *get_edge_indices(int edge, int ori, order1_t order) {
		CHECK_EDGE(edge);
		if (!edge_indices[edge][ori].exists(order)) compute_edge_indices(edge, ori, order);
		return edge_indices[edge][ori][order];
	}

	virtual int *get_face_indices(int face, int ori, order2_t order) {
		CHECK_FACE(face);
		if (!face_indices[face][ori].exists(order.get_idx())) compute_face_indices(face, ori, order);
		return face_indices[face][ori][order.get_idx()];
	}

  	virtual int *get_bubble_indices(order3_t order) {
		if (!bubble_indices.exists(order.get_idx())) compute_bubble_indices(order);
		return bubble_indices[order.get_idx()];
	}

	virtual int get_num_edge_fns(order1_t order) const {
		if (order > 1) return (order - 1);
		else return 0;
	}

	virtual int get_num_face_fns(order2_t order) const {
		if (order.x > 1 && order.y > 1) return (order.x - 1) * (order.y - 1);
		else return 0;
	}

	virtual int get_num_bubble_fns(order3_t order) const {
		if (order.x > 1 && order.y > 1 && order.z > 1) return (order.x - 1) * (order.y - 1) * (order.z - 1);
		else return 0;
	}

	virtual int get_face_orientations(int face) const { return RefHex::get_face_orientations(face); }

	virtual int get_edge_orientations() const { return RefHex::get_edge_orientations(); }

	virtual order3_t get_order(int index) const;

	virtual int get_shape_type(int index) const { return -1; }

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

	virtual double get_val(int n, int index, double x, double y, double z, int component) {
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

