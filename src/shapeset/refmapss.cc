// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
// - Pavel Kus
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
#include "../shapeset.h"
#include "common.h"
#include "refmapss.h"
#include <common/error.h>
#include <common/callstack.h>
#include "lobatto.h"

//// RefMapShapesetTetra //////////////////////////////////////////////////////////////////////////

#ifdef WITH_TETRA

// shape functions for tetrahedron

static double refmap_tetra_f0(double x, double y, double z) { return lambda1(x, y, z); }
static double refmap_tetra_f1(double x, double y, double z) { return lambda2(x, y, z); }
static double refmap_tetra_f2(double x, double y, double z) { return lambda0(x, y, z); }
static double refmap_tetra_f3(double x, double y, double z) { return lambda3(x, y, z); }

static shape_fn_t refmap_tetra_fn[] = {
	refmap_tetra_f0, refmap_tetra_f1, refmap_tetra_f2, refmap_tetra_f3
};

// DX /////////////////////////////////////////////////////////////////////////////////////////////

double refmap_tetra_dx_f0(double x, double y, double z) { return lambda1dx(x, y, z); }
double refmap_tetra_dx_f1(double x, double y, double z) { return lambda2dx(x, y, z); }
double refmap_tetra_dx_f2(double x, double y, double z) { return lambda0dx(x, y, z); }
double refmap_tetra_dx_f3(double x, double y, double z) { return lambda3dx(x, y, z); }

static shape_fn_t refmap_tetra_dx[] = {
	refmap_tetra_dx_f0, refmap_tetra_dx_f1, refmap_tetra_dx_f2, refmap_tetra_dx_f3
};

// DY /////////////////////////////////////////////////////////////////////////////////////////////

double refmap_tetra_dy_f0(double x, double y, double z) { return lambda1dy(x, y, z); }
double refmap_tetra_dy_f1(double x, double y, double z) { return lambda2dy(x, y, z); }
double refmap_tetra_dy_f2(double x, double y, double z) { return lambda0dy(x, y, z); }
double refmap_tetra_dy_f3(double x, double y, double z) { return lambda3dy(x, y, z); }

static shape_fn_t refmap_tetra_dy[] = {
	refmap_tetra_dy_f0, refmap_tetra_dy_f1, refmap_tetra_dy_f2, refmap_tetra_dy_f3
};

// DZ /////////////////////////////////////////////////////////////////////////////////////////////

double refmap_tetra_dz_f0(double x, double y, double z) { return lambda1dz(x, y, z); }
double refmap_tetra_dz_f1(double x, double y, double z) { return lambda2dz(x, y, z); }
double refmap_tetra_dz_f2(double x, double y, double z) { return lambda0dz(x, y, z); }
double refmap_tetra_dz_f3(double x, double y, double z) { return lambda3dz(x, y, z); }

static shape_fn_t refmap_tetra_dz[] = {
	refmap_tetra_dz_f0, refmap_tetra_dz_f1, refmap_tetra_dz_f2, refmap_tetra_dz_f3
};

//

static int refmap_tetra_vertex_indices[] = { 0, 1, 2, 3 };

static shape_fn_t *refmap_tetra_fn_table[] = { refmap_tetra_fn };
static shape_fn_t *refmap_tetra_dx_table[] = { refmap_tetra_dx };
static shape_fn_t *refmap_tetra_dy_table[] = { refmap_tetra_dy };
static shape_fn_t *refmap_tetra_dz_table[] = { refmap_tetra_dz };

#endif

RefMapShapesetTetra::RefMapShapesetTetra() : Shapeset(4)
{
#ifdef WITH_TETRA
	type = H1;
	mode = MODE_TETRAHEDRON;
	num_components = 1;

	// fn values are calculated by the tables
	shape_table[FN]  = refmap_tetra_fn_table;
	shape_table[DX]  = refmap_tetra_dx_table;
	shape_table[DY]  = refmap_tetra_dy_table;
	shape_table[DZ]  = refmap_tetra_dz_table;
	shape_table[DXY] = NULL;
	shape_table[DXZ] = NULL;
	shape_table[DYZ] = NULL;

	vertex_indices = refmap_tetra_vertex_indices;
#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}

RefMapShapesetTetra::~RefMapShapesetTetra() {
}

//// RefMapShapesetHex ////////////////////////////////////////////////////////////////////////////

#ifdef WITH_HEX

struct h1_hex_index_t {
	unsigned type:2;
	unsigned ef:4;
	unsigned ori:3;
	unsigned x:4;
	unsigned y:4;
	unsigned z:4;

	h1_hex_index_t(int idx) {
		this->type = (idx >> 19) & 0x03;
		this->ef =   (idx >> 15) & 0x0F;
		this->ori =  (idx >> 12) & 0x07;
		this->x = (idx >> 8) & 0x0F;
		this->y = (idx >> 4) & 0x0F;
		this->z = (idx >> 0) & 0x0F;
	}

	h1_hex_index_t(int type, int ef, int x, int y, int z, int ori = 0) {
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


static int lobatto_hex_vertex_indices[] = {
	h1_hex_index_t(SHFN_VERTEX, 0, 0, 0, 0), h1_hex_index_t(SHFN_VERTEX, 0, 1, 0, 0),
	h1_hex_index_t(SHFN_VERTEX, 0, 1, 1, 0), h1_hex_index_t(SHFN_VERTEX, 0, 0, 1, 0),
	h1_hex_index_t(SHFN_VERTEX, 0, 0, 0, 1), h1_hex_index_t(SHFN_VERTEX, 0, 1, 0, 1),
	h1_hex_index_t(SHFN_VERTEX, 0, 1, 1, 1), h1_hex_index_t(SHFN_VERTEX, 0, 0, 1, 1)
};

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


static void decompose(h1_hex_index_t index, int indices[3], int ori[3], bool swapori = true) {
	_F_
	int permut[3];
	int num_01;

//	printf("idx = %d\n", (int) index);

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

static void calc_fn_values(int index, int np, QuadPt3D *pt, int component, double *val) {
	_F_
	h1_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int k = 0; k < np; k++) {
		double x = (oris[0] == 0) ? pt[k].x : -pt[k].x;
		double y = (oris[1] == 0) ? pt[k].y : -pt[k].y;
		double z = (oris[2] == 0) ? pt[k].z : -pt[k].z;

		val[k] = lobatto_fn_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
	}
}


static void calc_dx_values(int index, int np, QuadPt3D *pt, int component, double *dx) {
	_F_
	h1_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	for (int k = 0; k < np; k++) {
		double x = (oris[0] == 0) ? pt[k].x : -pt[k].x;
		double y = (oris[1] == 0) ? pt[k].y : -pt[k].y;
		double z = (oris[2] == 0) ? pt[k].z : -pt[k].z;

		dx[k] = lobatto_der_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
		if (oris[0] == 1) dx[k] = -dx[k];
	}
}


static void calc_dy_values(int index, int np, QuadPt3D *pt, int component, double *dy) {
	_F_
	h1_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	for (int k = 0; k < np; k++) {
		double x = (oris[0] == 0) ? pt[k].x : -pt[k].x;
		double y = (oris[1] == 0) ? pt[k].y : -pt[k].y;
		double z = (oris[2] == 0) ? pt[k].z : -pt[k].z;

		dy[k] = lobatto_fn_tab_1d[indices[0]](x) * lobatto_der_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
		if (oris[1] == 1) dy[k] = -dy[k];
	}
}


static void calc_dz_values(int index, int np, QuadPt3D *pt, int component, double *dz) {
	_F_
	h1_hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	for (int k = 0; k < np; k++) {
		double x = (oris[0] == 0) ? pt[k].x : -pt[k].x;
		double y = (oris[1] == 0) ? pt[k].y : -pt[k].y;
		double z = (oris[2] == 0) ? pt[k].z : -pt[k].z;

		dz[k] = lobatto_fn_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_der_tab_1d[indices[2]](z);
		if (oris[2] == 1) dz[k] = -dz[k];
	}
}

/*
// shape functions for hexahedron

static double refmap_hex_f0(double x, double y, double z) { return l0(x) * l0(y) * l0(z); }
static double refmap_hex_f1(double x, double y, double z) { return l1(x) * l0(y) * l0(z); }
static double refmap_hex_f2(double x, double y, double z) { return l1(x) * l1(y) * l0(z); }
static double refmap_hex_f3(double x, double y, double z) { return l0(x) * l1(y) * l0(z); }
static double refmap_hex_f4(double x, double y, double z) { return l0(x) * l0(y) * l1(z); }
static double refmap_hex_f5(double x, double y, double z) { return l1(x) * l0(y) * l1(z); }
static double refmap_hex_f6(double x, double y, double z) { return l1(x) * l1(y) * l1(z); }
static double refmap_hex_f7(double x, double y, double z) { return l0(x) * l1(y) * l1(z); }

static shape_fn_t refmap_hex_fn[] = {
	refmap_hex_f0, refmap_hex_f1, refmap_hex_f2, refmap_hex_f3,
	refmap_hex_f4, refmap_hex_f5, refmap_hex_f6, refmap_hex_f7
};

// DX /////////////////////////////////////////////////////////////////////////////////////////////

static double refmap_hex_dx_f0(double x, double y, double z) { return dl0(x) * l0(y) * l0(z); }
static double refmap_hex_dx_f1(double x, double y, double z) { return dl1(x) * l0(y) * l0(z); }
static double refmap_hex_dx_f2(double x, double y, double z) { return dl1(x) * l1(y) * l0(z); }
static double refmap_hex_dx_f3(double x, double y, double z) { return dl0(x) * l1(y) * l0(z); }
static double refmap_hex_dx_f4(double x, double y, double z) { return dl0(x) * l0(y) * l1(z); }
static double refmap_hex_dx_f5(double x, double y, double z) { return dl1(x) * l0(y) * l1(z); }
static double refmap_hex_dx_f6(double x, double y, double z) { return dl1(x) * l1(y) * l1(z); }
static double refmap_hex_dx_f7(double x, double y, double z) { return dl0(x) * l1(y) * l1(z); }

static shape_fn_t refmap_hex_dx[] = {
	refmap_hex_dx_f0, refmap_hex_dx_f1, refmap_hex_dx_f2, refmap_hex_dx_f3,
	refmap_hex_dx_f4, refmap_hex_dx_f5, refmap_hex_dx_f6, refmap_hex_dx_f7
};

// DY /////////////////////////////////////////////////////////////////////////////////////////////

static double refmap_hex_dy_f0(double x, double y, double z) { return l0(x) * dl0(y) * l0(z); }
static double refmap_hex_dy_f1(double x, double y, double z) { return l1(x) * dl0(y) * l0(z); }
static double refmap_hex_dy_f2(double x, double y, double z) { return l1(x) * dl1(y) * l0(z); }
static double refmap_hex_dy_f3(double x, double y, double z) { return l0(x) * dl1(y) * l0(z); }
static double refmap_hex_dy_f4(double x, double y, double z) { return l0(x) * dl0(y) * l1(z); }
static double refmap_hex_dy_f5(double x, double y, double z) { return l1(x) * dl0(y) * l1(z); }
static double refmap_hex_dy_f6(double x, double y, double z) { return l1(x) * dl1(y) * l1(z); }
static double refmap_hex_dy_f7(double x, double y, double z) { return l0(x) * dl1(y) * l1(z); }

static shape_fn_t refmap_hex_dy[] = {
	refmap_hex_dy_f0, refmap_hex_dy_f1, refmap_hex_dy_f2, refmap_hex_dy_f3,
	refmap_hex_dy_f4, refmap_hex_dy_f5, refmap_hex_dy_f6, refmap_hex_dy_f7
};

// DZ /////////////////////////////////////////////////////////////////////////////////////////////

static double refmap_hex_dz_f0(double x, double y, double z) { return l0(x) * l0(y) * dl0(z); }
static double refmap_hex_dz_f1(double x, double y, double z) { return l1(x) * l0(y) * dl0(z); }
static double refmap_hex_dz_f2(double x, double y, double z) { return l1(x) * l1(y) * dl0(z); }
static double refmap_hex_dz_f3(double x, double y, double z) { return l0(x) * l1(y) * dl0(z); }
static double refmap_hex_dz_f4(double x, double y, double z) { return l0(x) * l0(y) * dl1(z); }
static double refmap_hex_dz_f5(double x, double y, double z) { return l1(x) * l0(y) * dl1(z); }
static double refmap_hex_dz_f6(double x, double y, double z) { return l1(x) * l1(y) * dl1(z); }
static double refmap_hex_dz_f7(double x, double y, double z) { return l0(x) * l1(y) * dl1(z); }

static shape_fn_t refmap_hex_dz[] = {
	refmap_hex_dz_f0, refmap_hex_dz_f1, refmap_hex_dz_f2, refmap_hex_dz_f3,
	refmap_hex_dz_f4, refmap_hex_dz_f5, refmap_hex_dz_f6, refmap_hex_dz_f7
};

//

static int refmap_hex_vertex_indices[] = { 0, 1, 2, 3, 4, 5, 6, 7 };

static shape_fn_t *refmap_hex_fn_table[] = { refmap_hex_fn };
static shape_fn_t *refmap_hex_dx_table[] = { refmap_hex_dx };
static shape_fn_t *refmap_hex_dy_table[] = { refmap_hex_dy };
static shape_fn_t *refmap_hex_dz_table[] = { refmap_hex_dz };
*/

#endif


RefMapShapesetHex::RefMapShapesetHex() : Shapeset(5)
{
	_F_
#ifdef WITH_HEX
/*	type = H1;
	mode = MODE_HEXAHEDRON;
	num_components = 1;

	// nothing is tabelated
	shape_table[FN]  = refmap_hex_fn_table;
	shape_table[DX]  = refmap_hex_dx_table;
	shape_table[DY]  = refmap_hex_dy_table;
	shape_table[DZ]  = refmap_hex_dz_table;
	shape_table[DXY] = NULL;
	shape_table[DXZ] = NULL;
	shape_table[DYZ] = NULL;

	// vertices
	vertex_indices = refmap_hex_vertex_indices;
*/
	type = H1;
	mode = MODE_HEXAHEDRON;
	num_components = 1;

	// fn, dx, dy, dz will be calculated on-the-fly
	shape_table_deleg[FN]  = calc_fn_values;
	shape_table_deleg[DX]  = calc_dx_values;
	shape_table_deleg[DY]  = calc_dy_values;
	shape_table_deleg[DZ]  = calc_dz_values;
	shape_table_deleg[DXY] = NULL;
	shape_table_deleg[DXZ] = NULL;
	shape_table_deleg[DYZ] = NULL;

	// vertices
	vertex_indices = lobatto_hex_vertex_indices;

#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}

RefMapShapesetHex::~RefMapShapesetHex() {
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

order3_t RefMapShapesetHex::get_order(int index) const {
	_F_
#ifdef WITH_HEX
	if (index >= 0) {
		h1_hex_index_t idx(index);
		order3_t ord(lobatto_order_1d[idx.x], lobatto_order_1d[idx.y], lobatto_order_1d[idx.z]);
		if (idx.type == SHFN_FACE && idx.ori >= 4) ord = turn_hex_face_order(idx.ef, ord);		// face function is turned due to orientation
		return ord;
	}
	else
		return get_ced_order(index);
#else
	return order3_t(0, 0, 0);
#endif
}

order3_t RefMapShapesetHex::get_dcmp(int index) const
{
	if (index >= 0) {
		h1_hex_index_t idx(index);
		order3_t ord(idx.x, idx.y, idx.z);
		return ord;
	}
	else
		return order3_t(-1);
}

int RefMapShapesetHex::get_shape_type(int index) const
{
	_F_
#ifdef WITH_HEX
	if (index >= 0) {
		h1_hex_index_t idx(index);
		return idx.type;
	}
	else
		return SHFN_CONSTRAINED;
#else
	return -1;
#endif
}

void RefMapShapesetHex::compute_edge_indices(int edge, int ori, order1_t order) {
	_F_
#ifdef WITH_HEX
	assert(order > 1);
	int *indices = new int[order - 1];
	MEM_CHECK(indices);

	int idx = 0;
	switch (edge) {
		case  0: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 0, i, 0, 0, ori); break;
		case  1: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 1, 1, i, 0, ori); break;
		case  2: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 2, i, 1, 0, ori); break;
		case  3: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 3, 0, i, 0, ori); break;
		case  4: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 4, 0, 0, i, ori); break;
		case  5: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 5, 1, 0, i, ori); break;
		case  6: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 6, 1, 1, i, ori); break;
		case  7: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 7, 0, 1, i, ori); break;
		case  8: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 8, i, 0, 1, ori); break;
		case  9: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 9, 1, i, 1, ori); break;
		case 10: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 10, i, 1, 1, ori); break;
		case 11: for (int i = 2; i <= order; i++) indices[idx++] = h1_hex_index_t(SHFN_EDGE, 11, 0, i, 1, ori); break;
		default: EXIT("Invalid edge number %d. Can be 0 - 11.", edge); break;
	}

	edge_indices[edge][ori][order] = indices;
#endif
}

void RefMapShapesetHex::compute_face_indices(int face, int ori, order2_t order) {
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
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 0, 0, i, j, ori);
			break;

		case 1:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 1, 1, i, j, ori);
			break;

		case 2:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 2, i, 0, j, ori);
			break;

		case 3:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 3, i, 1, j, ori);
			break;

		case 4:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 4, i, j, 0, ori);
			break;

		case 5:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = h1_hex_index_t(SHFN_FACE, 5, i, j, 1, ori);
			break;

		default:
			EXIT("Invalid face number %d. Can be 0 - 5.", face);
			break;
	}

	face_indices[face][ori][order.get_idx()] = indices;
#endif
}

void RefMapShapesetHex::compute_bubble_indices(order3_t order) {
	_F_
#ifdef WITH_HEX
	assert(order.x > 1);
	assert(order.y > 1);
	assert(order.z > 1);
	int *indices = new int[(order.x - 1) * (order.y - 1) * (order.z - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	for (unsigned int i = 2; i <= order.x; i++)
		for (unsigned int j = 2; j <= order.y; j++)
			for (unsigned int k = 2; k <= order.z; k++) {
				indices[idx++] = h1_hex_index_t(SHFN_BUBBLE, 0, i, j, k, 0);
//				printf("%d\n", indices[idx - 1]);
			}

	bubble_indices[order.get_idx()] = indices;
#endif
}
