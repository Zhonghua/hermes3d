// This file is part of Hermes3D
//
// Copyright (c) 2007 - 2009 David Andrs <dandrs@unr.edu>
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

//
// h1sinhex.cc
//

#include "../h3dconfig.h"
#include "common.h"
#include "h1sinhex.h"
#include <common/error.h>
#include <common/callstack.h>
#include "matrix.h"
#include "lobatto.h"

#include "mesh.h"

#ifdef WITH_HEX

struct h1s_hex_index_t {
	unsigned type:2;
	unsigned ef:4;
	unsigned ori:3;
	unsigned x:4;
	unsigned y:4;
	unsigned z:4;

	h1s_hex_index_t(int idx) {
		this->type = (idx >> 19) & 0x03;
		this->ef =   (idx >> 15) & 0x0F;
		this->ori =  (idx >> 12) & 0x07;
		this->x = (idx >> 8) & 0x0F;
		this->y = (idx >> 4) & 0x0F;
		this->z = (idx >> 0) & 0x0F;
	}

	h1s_hex_index_t(int type, int ef, int x, int y, int z, int ori = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->ori = ori;
		this->type = type;
		this->ef = ef;
	}

	operator int() { return (type << 19) | (ef << 15) | (ori << 12) | (x << 8) | (y << 4) | z; }

	void to_str() { printf("hex = [type=%d, ef=%d, ori=%d, x=%d, y=%d, z=%d]", type, ef, ori, x, y, z); }
};

static scalar fn_0(double x)  { return l0(x); }
static scalar fn_1(double x)  { return l1(x); }
static scalar fn_2(double x)  { return l2(x); }
static scalar fn_3(double x)  { return l3(x); }
static scalar fn_4(double x)  { return l4(x); }
static scalar fn_5(double x)  { return l5(x); }
static scalar fn_6(double x)  { return l6(x); }
static scalar fn_7(double x)  { return l7(x); }
static scalar fn_8(double x)  { return l8(x); }
static scalar fn_9(double x)  { return l9(x); }
static scalar fn_10(double x) { return l10(x); }

static shape_fn_1d_t fn_tab_1d[] = { fn_0, fn_1, fn_2, fn_3, fn_4, fn_5, fn_6, fn_7, fn_8, fn_9, fn_10 };

static scalar der_0(double x)  { return dl0(x); }
static scalar der_1(double x)  { return dl1(x); }
static scalar der_2(double x)  { return dl2(x); }
static scalar der_3(double x)  { return dl3(x); }
static scalar der_4(double x)  { return dl4(x); }
static scalar der_5(double x)  { return dl5(x); }
static scalar der_6(double x)  { return dl6(x); }
static scalar der_7(double x)  { return dl7(x); }
static scalar der_8(double x)  { return dl8(x); }
static scalar der_9(double x)  { return dl9(x); }
static scalar der_10(double x) { return dl10(x); }

static shape_fn_1d_t der_tab_1d[] = { der_0, der_1, der_2, der_3, der_4, der_5, der_6, der_7, der_8, der_9, der_10 };

static int hex_vertex_indices[] = {
	h1s_hex_index_t(SHFN_VERTEX, 0, 0, 0, 0), h1s_hex_index_t(SHFN_VERTEX, 0, 1, 0, 0),
	h1s_hex_index_t(SHFN_VERTEX, 0, 1, 1, 0), h1s_hex_index_t(SHFN_VERTEX, 0, 0, 1, 0),
	h1s_hex_index_t(SHFN_VERTEX, 0, 0, 0, 1), h1s_hex_index_t(SHFN_VERTEX, 0, 1, 0, 1),
	h1s_hex_index_t(SHFN_VERTEX, 0, 1, 1, 1), h1s_hex_index_t(SHFN_VERTEX, 0, 0, 1, 1)
};

static int index_order[] = { 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

// -- helpers -- //

static void find_permutation(int *indices, int *permut, int &num_01) {
	_F_
	for (int i = 0; i < 3; i++)
		permut[i] = i;
	num_01 = 0;
	if (indices[0] < 2) num_01++;
	if (indices[1] < 2) {
		num_01++;
		if (num_01 == 1) std::swap(permut[0], permut[1]);
	}
	if (indices[2] < 2) {
		num_01++;
		if (num_01 == 1) {
			std::swap(permut[1], permut[2]);
			std::swap(permut[0], permut[1]);
		}
		if (num_01 == 2) std::swap(permut[1], permut[2]);
	}
}


static void decompose(h1s_hex_index_t index, int indices[3], int ori[3], bool swapori = true) {
	_F_
	int permut[3];
	int num_01;

	// get order in each direction
	indices[0] = index.x;
	indices[1] = index.y;
	indices[2] = index.z;

	find_permutation(indices, permut, num_01);

	ori[0] = ori[1] = ori[2] = 0;
	if (num_01 == 2) {
		// edge function
		assert(index.ori == 0 || index.ori == 1);
		ori[permut[2]] = index.ori;
	}
	else if (num_01 == 1) {
		// face function
		assert(index.ori >= 0 && index.ori < 8);
		if (index.ori % 2 == 1) ori[permut[1]] = 1;
		if (index.ori % 4 >= 2) ori[permut[2]] = 1;
		if (index.ori >= 4) {
			std::swap(indices[permut[1]], indices[permut[2]]);
			if (swapori) std::swap(ori[permut[1]], ori[permut[2]]);
		}
	}
	else {
		// vertex or bubble function
		assert(index.ori == 0);
	}
}

// -- functions that calculate values of fn, dx, dy, dz on the fly -- //

static double calc_fn_value(int index, double x, double y, double z, int component) {
	_F_
	h1s_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	return lobatto_fn_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
}


static double calc_dx_value(int index, double x, double y, double z, int component) {
	_F_
	h1s_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	double dx = lobatto_der_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
	if (oris[0] == 1) dx = -dx;

	return dx;
}


static double calc_dy_value(int index, double x, double y, double z, int component) {
	_F_
	h1s_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	double dy = lobatto_fn_tab_1d[indices[0]](x) * lobatto_der_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
	if (oris[1] == 1) dy = -dy;

	return dy;
}


static double calc_dz_value(int index, double x, double y, double z, int component) {
	_F_
	h1s_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	double dz = lobatto_fn_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_der_tab_1d[indices[2]](z);
	if (oris[2] == 1) dz = -dz;

	return dz;
}

#endif

// -- -- //

H1ShapesetSinHex::H1ShapesetSinHex() {
	_F_
#ifdef WITH_HEX
	type = H1;
	mode = MODE_HEXAHEDRON;
	num_components = 1;

	// fn, dx, dy, dz will be calculated on-the-fly
	shape_table_deleg[FN]  = calc_fn_value;
	shape_table_deleg[DX]  = calc_dx_value;
	shape_table_deleg[DY]  = calc_dy_value;
	shape_table_deleg[DZ]  = calc_dz_value;
	shape_table_deleg[DXY] = NULL;
	shape_table_deleg[DXZ] = NULL;
	shape_table_deleg[DYZ] = NULL;

	// vertices
	vertex_indices = hex_vertex_indices;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

H1ShapesetSinHex::~H1ShapesetSinHex() {
	_F_
#ifdef WITH_HEX
	for (int edge = 0; edge < Hex::NUM_EDGES; edge++)
		for (int ori = 0; ori < NUM_EDGE_ORIS; ori++)
			for (Word_t idx = edge_indices[edge][ori].first(); idx != INVALID_IDX; idx = edge_indices[edge][ori].next(idx))
				delete [] edge_indices[edge][ori][idx];

	for (int face = 0; face < Hex::NUM_FACES; face++)
		for (int ori = 0; ori < NUM_FACE_ORIS; ori++)
			for (Word_t idx = face_indices[face][ori].first(); idx != INVALID_IDX; idx = face_indices[face][ori].next(idx))
				delete [] face_indices[face][ori][idx];

	for (Word_t idx = bubble_indices.first(); idx != INVALID_IDX; idx = bubble_indices.next(idx))
		delete [] bubble_indices[idx];
#endif
}

order3_t H1ShapesetSinHex::get_order(int index) const {
	_F_
#ifdef WITH_HEX
	if (index >= 0) {
		h1s_hex_index_t idx(index);
		order3_t ord(index_order[idx.x], index_order[idx.y], index_order[idx.z]);
		if (idx.ori >= 4) ord = turn_hex_face_order(idx.ef, ord);		// face function is turned due to orientation
		return ord;
	}
	else {
		EXIT(ERR_NOT_IMPLEMENTED);
		return order3_t(0, 0, 0);
	}
#else
	return order3_t(0, 0, 0);
#endif
}

void H1ShapesetSinHex::compute_edge_indices(int edge, int ori, order1_t order) {
	_F_
#ifdef WITH_HEX
	assert(order > 1);
	int *indices = new int[order - 1];
	MEM_CHECK(indices);

	int idx = 0;
	switch (edge) {
		case  0: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 0, i, 0, 0, ori); break;
		case  1: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 1, 1, i, 0, ori); break;
		case  2: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 2, i, 1, 0, ori); break;
		case  3: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 3, 0, i, 0, ori); break;
		case  4: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 4, 0, 0, i, ori); break;
		case  5: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 5, 1, 0, i, ori); break;
		case  6: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 6, 1, 1, i, ori); break;
		case  7: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 7, 0, 1, i, ori); break;
		case  8: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 8, i, 0, 1, ori); break;
		case  9: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 9, 1, i, 1, ori); break;
		case 10: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 10, i, 1, 1, ori); break;
		case 11: for (int i = 2; i <= order; i++) indices[idx++] = h1s_hex_index_t(SHFN_EDGE, 11, 0, i, 1, ori); break;
		default: EXIT(ERR_FAILURE, "Invalid edge number %d. Can be 0 - 11.", edge); break;
	}

	edge_indices[edge][ori][order] = indices;
#endif
}

void H1ShapesetSinHex::compute_face_indices(int face, int ori, order2_t order) {
	_F_
#ifdef WITH_HEX
	assert(order.x > 1);
	assert(order.y > 1);
	int horder = order.x, vorder = order.y;
	int *indices = new int[(horder - 1) * (vorder - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	switch (face) {
		case 0:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1s_hex_index_t(SHFN_FACE, 0, 0, i, j, ori);
			break;

		case 1:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1s_hex_index_t(SHFN_FACE, 1, 1, i, j, ori);
			break;

		case 2:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1s_hex_index_t(SHFN_FACE, 2, i, 0, j, ori);
			break;

		case 3:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1s_hex_index_t(SHFN_FACE, 3, i, 1, j, ori);
			break;

		case 4:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1s_hex_index_t(SHFN_FACE, 4, i, j, 0, ori);
			break;

		case 5:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1s_hex_index_t(SHFN_FACE, 5, i, j, 1, ori);
			break;

		default:
			EXIT(ERR_FAILURE, "Invalid face number %d. Can be 0 - 5.", face);
			break;
	}

	face_indices[face][ori][order.get_idx()] = indices;
#endif
}

void H1ShapesetSinHex::compute_bubble_indices(order3_t order) {
	_F_
#ifdef WITH_HEX
	assert(order.x > 1);
	assert(order.y > 1);
	assert(order.z > 1);
	int *indices = new int[(order.x - 1) * (order.y - 1) * (order.z - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	for (int i = 2; i <= order.x; i++)
		for(int j = 2; j <= order.y; j++)
			for(int k = 2; k <= order.z; k++)
				indices[idx++] = h1s_hex_index_t(SHFN_BUBBLE, 0, i, j, k, 0);

	bubble_indices[order.get_idx()] = indices;
#endif
}

/// --- CED specific stuff ---

CEDComb *H1ShapesetSinHex::calc_constrained_edge_combination(int ori, int order, Part part) {
	_F_
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}

CEDComb *H1ShapesetSinHex::calc_constrained_edge_face_combination(int ori, int order, Part part, int dir) {
	_F_
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}

CEDComb *H1ShapesetSinHex::calc_constrained_face_combination(int ori, int order, Part part) {
	_F_
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}
