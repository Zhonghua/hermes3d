// This file is part of Hermes3D
//
// Copyright (c) 2008 - 2009 Pavel Kus <pavel.kus@gmail.com>
// Copyright (c) 2009 David Andrs <dandrs@unr.edu>
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

#include "../h3dconfig.h"
#include "../refdomain.h"
#include "common.h"
#include "lobatto.h"
#include "hcurllobattohex.h"
#include <common/error.h>
#include <common/callstack.h>
#include "matrix.h"

#include "mesh.h"

#ifdef WITH_HEX

struct hc_hex_index_t {
	unsigned type:2;		// SHFN_XXX
	unsigned ef:4;			// edge/face index (depends on type)
	unsigned ori:3;
	unsigned l:2;			// legendre (0-2)
	unsigned x:4;
	unsigned y:4;
	unsigned z:4;

	hc_hex_index_t(int idx) {
		this->x = (idx >> 8) & 0x0F;
		this->y = (idx >> 4) & 0x0F;
		this->z = (idx >> 0) & 0x0F;

		this->l    = (idx >> 12) & 0x03;
		this->ori  = (idx >> 14) & 0x07;
		this->ef   = (idx >> 17) & 0x0F;
		this->type = (idx >> 21) & 0x03;
	}

	hc_hex_index_t(int type, int ef, int x, int y, int z, int l, int ori = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->l = l;
		this->ori = ori;
		this->ef = ef;
		this->type = type;
	}

	operator int() { return (type << 21) | (ef << 17) | (ori << 14) | (l << 12) | (x  << 8) | (y << 4) | z; }
};


static void decompose(hc_hex_index_t ind, int indices[3], int ori[3], int &which_legendre) {
	_F_
	indices[0] = ind.x;
	indices[1] = ind.y;
	indices[2] = ind.z;

	ori[0] = 0;
	ori[1] = 0;
	ori[2] = 0;

	which_legendre = ind.l;

	if (ind.type == SHFN_EDGE) {
		assert(ind.ori == 0 || ind.ori == 1);
		ori[ind.l] = ind.ori;
	}
	else if (ind.type == SHFN_FACE) {
		assert(ind.ori >= 0 && ind.ori <= 7);
		int dir_1 = RefHex::get_face_tangent_direction(ind.ef, 0);
		int dir_2 = RefHex::get_face_tangent_direction(ind.ef, 1);

		if (ind.ori % 2 == 1) ori[dir_1] = 1;
		if (ind.ori % 4 >= 2) ori[dir_2] = 1;
		if (ind.ori >= 4) {
			std::swap(indices[dir_1], indices[dir_2]);
			std::swap(ori[dir_1], ori[dir_2]);
			which_legendre = (which_legendre == dir_1) ? dir_2 : dir_1;
		}
	}
	else {
		assert(ind.ori == 0);
	}
}

static void calc_fn_values(int index, int np, QuadPt3D *pt, int component, double *val) {
	_F_
	hc_hex_index_t ind(index);
	int indices[3];
	int ori[3];
	int which_legendre;

	decompose(ind, indices, ori, which_legendre);

	if (which_legendre == component) {
		for (int k = 0; k < np; k++) {
			double point[3] = { pt[k].x, pt[k].y, pt[k].z };

			if (ori[0] == 1) point[0] = -point[0];
			if (ori[1] == 1) point[1] = -point[1];
			if (ori[2] == 1) point[2] = -point[2];

			double value = 1.0;
			for (int i = 0; i < 3; i++) {
				if (which_legendre == i) value *= legendre_fn_tab_1d[indices[i]](point[i]);
				else value *= lobatto_fn_tab_1d[indices[i]](point[i]);
			}

			if (ori[component] == 1) value = -value;
			val[k] = value;
		}
	}
	else {
		for (int k = 0; k < np; k++)
			val[k] = 0.0;
	}
}

static void calc_der_values(int index, int np, QuadPt3D *pt, int component, int which_der, double *val) {
	_F_
	hc_hex_index_t ind(index);
	int indices[3];
	int ori[3];
	int which_legendre;

	decompose(ind, indices, ori, which_legendre);

	// for the sake of simplicity, recode orientations like this: 0 -> 1, 1 -> -1
	for (int i = 0; i < 3; i++)
		ori[i] = (ori[i] == 0) ? 1 : -1;

	if (which_legendre == component) {
		for (int k = 0; k < np; k++) {
			double point[3] = { pt[k].x, pt[k].y, pt[k].z };

			if (ori[0] == -1) point[0] = -point[0];
			if (ori[1] == -1) point[1] = -point[1];
			if (ori[2] == -1) point[2] = -point[2];

			double value = 1.0;
			for (int i = 0; i < 3; i++) {
				if (which_legendre == i) {
					if (which_der == i) value *= legendre_der_tab_1d[indices[i]](point[i]);
					else value *= legendre_fn_tab_1d[indices[i]](point[i]);
				}
				else {
					if (which_der == i) value *= lobatto_der_tab_1d[indices[i]](point[i]);
					else value *= lobatto_fn_tab_1d[indices[i]](point[i]);
				}
			}

			val[k] = ori[component] * ori[which_der] * value;
		}
	}
	else {
		for (int k = 0; k < np; k++) val[k] = 0.0;
	}
}

static void calc_dx_values(int index, int np, QuadPt3D *pt, int component, double *val) {
	_F_
	return calc_der_values(index, np, pt, component, 0, val);
}

static void calc_dy_values(int index, int np, QuadPt3D *pt, int component, double *val) {
	_F_
	return calc_der_values(index, np, pt, component, 1, val);
}

static void calc_dz_values(int index, int np, QuadPt3D *pt, int component, double *val) {
	_F_
	return calc_der_values(index, np, pt, component, 2, val);
}

#endif

HcurlShapesetLobattoHex::HcurlShapesetLobattoHex() {
	_F_
#ifdef WITH_HEX
	type = Hcurl;
	mode = MODE_HEXAHEDRON;

	num_components = 3;

	// fn, dx, dy, dz will be calculated on-the-fly
	shape_table_deleg[FN]  = calc_fn_values;
	shape_table_deleg[DX]  = calc_dx_values;
	shape_table_deleg[DY]  = calc_dy_values;
	shape_table_deleg[DZ]  = calc_dz_values;
	shape_table_deleg[DXY] = NULL;
	shape_table_deleg[DXZ] = NULL;
	shape_table_deleg[DYZ] = NULL;

	vertex_indices = NULL;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

HcurlShapesetLobattoHex::~HcurlShapesetLobattoHex() {
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

order3_t HcurlShapesetLobattoHex::get_order(int index) const {
	_F_
#ifdef WITH_HEX
	if (index >= 0) {
		order3_t ord;
		hc_hex_index_t idx(index);
		if (idx.l == 0) ord = order3_t(legendre_order_1d[idx.x], lobatto_order_1d[idx.y], lobatto_order_1d[idx.z]);
		else if (idx.l == 1) ord = order3_t(lobatto_order_1d[idx.x], legendre_order_1d[idx.y], lobatto_order_1d[idx.z]);
		else if (idx.l == 2) ord = order3_t(lobatto_order_1d[idx.x], lobatto_order_1d[idx.y], legendre_order_1d[idx.z]);

		if (idx.type == SHFN_FACE && idx.ori >= 4) ord = turn_hex_face_order(idx.ef, ord);		// face function is turned due to orientation
		return ord;
	}
	else {
		EXIT(ERR_NOT_IMPLEMENTED);
		return order3_t(-1);
	}
#endif
}

void HcurlShapesetLobattoHex::compute_edge_indices(int edge, int ori, order1_t order) {
	_F_
#ifdef WITH_HEX
	int *indices = new int[get_num_edge_fns(order)];
	MEM_CHECK(indices);

	int idx = 0;
	switch (edge) {
		case  0: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, i, 0, 0, 0, ori); break;
		case  1: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, 1, i, 0, 1, ori); break;
		case  2: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, i, 1, 0, 0, ori); break;
		case  3: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, 0, i, 0, 1, ori); break;
		case  4: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, 0, 0, i, 2, ori); break;
		case  5: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, 1, 0, i, 2, ori); break;
		case  6: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, 1, 1, i, 2, ori); break;
		case  7: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, 0, 1, i, 2, ori); break;
		case  8: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, i, 0, 1, 0, ori); break;
		case  9: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, 1, i, 1, 1, ori); break;
		case 10: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, i, 1, 1, 0, ori); break;
		case 11: for (int i = 0; i <= order; i++) indices[idx++] = hc_hex_index_t(SHFN_EDGE, edge, 0, i, 1, 1, ori); break;
		default: EXIT(ERR_FAILURE, "Invalid edge number %d. Can be 0 - 11.", edge); break;
	}

	edge_indices[edge][ori][order] = indices;
#endif
}

void HcurlShapesetLobattoHex::compute_face_indices(int face, int ori, order2_t order) {
	_F_
#ifdef WITH_HEX
	int *indices = new int[get_num_face_fns(order)];
	MEM_CHECK(indices);

	int idx = 0;
	switch (face) {
		case 0:
			for (int i = 0; i <= order.x; i++)
				for (int j = 2; j <= order.y + 1; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, 0, i, j, 1, ori);
			for (int i = 2; i <= order.x + 1; i++)
				for (int j = 0; j <= order.y; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, 0, i, j, 2, ori);
			break;

		case 1:
			for (int i = 0; i <= order.x; i++)
				for (int j = 2; j <= order.y + 1; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, 1, i, j, 1, ori);
			for (int i = 2; i <= order.x + 1; i++)
				for (int j = 0; j <= order.y; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, 1, i, j, 2, ori);
			break;

		case 2:
			for (int i = 0; i <= order.x; i++)
				for (int j = 2; j <= order.y + 1; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, i, 0, j, 0, ori);
			for (int i = 2; i <= order.x + 1; i++)
				for (int j = 0; j <= order.y; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, i, 0, j, 2, ori);
			break;

		case 3:
			for (int i = 0; i <= order.x; i++)
				for (int j = 2; j <= order.y + 1; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, i, 1, j, 0, ori);
			for (int i = 2; i <= order.x + 1; i++)
				for (int j = 0; j <= order.y; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, i, 1, j, 2, ori);
			break;

		case 4:
			for (int i = 0; i <= order.x; i++)
				for (int j = 2; j <= order.y + 1; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, i, j, 0, 0, ori);
			for (int i = 2; i <= order.x + 1; i++)
				for (int j = 0; j <= order.y; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, i, j, 0, 1, ori);
			break;

		case 5:
			for (int i = 0; i <= order.x; i++)
				for (int j = 2; j <= order.y + 1; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, i, j, 1, 0, ori);
			for (int i = 2; i <= order.x + 1; i++)
				for (int j = 0; j <= order.y; j++)
					indices[idx++] = hc_hex_index_t(SHFN_FACE, face, i, j, 1, 1, ori);
			break;

		default:
			EXIT(ERR_FAILURE, "Invalid face number %d. Can be 0 - 5.", face);
			break;
	}

	face_indices[face][ori][order.get_idx()] = indices;
#endif
}

void HcurlShapesetLobattoHex::compute_bubble_indices(order3_t order) {
	_F_
#ifdef WITH_HEX
	int *indices = new int[get_num_bubble_fns(order)];
	MEM_CHECK(indices);

	int idx = 0;
	for (int i = 0; i <= order.x; i++)
		for (int j = 2; j <= order.y + 1; j++)
			for (int k = 2; k <= order.z + 1; k++)
				indices[idx++] = hc_hex_index_t(SHFN_BUBBLE, 0, i, j, k, 0);
	for (int i = 2; i <= order.x + 1; i++)
		for (int j = 0; j <= order.y; j++)
			for (int k = 2; k <= order.z + 1; k++)
				indices[idx++] = hc_hex_index_t(SHFN_BUBBLE, 0, i, j, k, 1);
	for (int i = 2; i <= order.x + 1; i++)
		for (int j = 2; j <= order.y + 1; j++)
			for (int k = 0; k <= order.z; k++)
				indices[idx++] = hc_hex_index_t(SHFN_BUBBLE, 0, i, j, k, 2);

	bubble_indices[order.get_idx()] = indices;
#endif
}

CEDComb *HcurlShapesetLobattoHex::calc_constrained_edge_combination(int ori, order1_t order, Part part) {
	_F_
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}

CEDComb *HcurlShapesetLobattoHex::calc_constrained_edge_face_combination(int ori, order2_t order, Part part) {
	_F_
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}

CEDComb *HcurlShapesetLobattoHex::calc_constrained_face_combination(int ori, order2_t order, Part part) {
	_F_
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}

