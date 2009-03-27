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

#include "config.h"
#include "shapeset.h"
#include "refdomain.h"
#include <common/trace.h>
#include <common/error.h>

/// TODO: move to common/shapeset?

int combine_face_part(int part, int finer_part) {
	assert(finer_part == 0 || finer_part == 1 || finer_part == 2);

	if (finer_part == 0) return part;							// stay the same
	else if (finer_part == 1) return get_lower_part(part);		// lower part
	else return get_higher_part(part);							// upper part
}

int opposite_part(int part) {
	int n;
	int m = part;
	for (n = 1; n <= part; n <<= 1)
		part -= n;
	return 3 * (n - 1) - m;
}

/// Get the endpoints of the interval
/// @param part[in] part of the interval (the number of stripe, see PIC)
/// @param lo[out] lower bound of the part the interval
/// @param hi[out] higher bound of the part the interval
void get_interval_part(int part, double &lo, double &hi) {
	int n;											// number of pieces of the interval
	for (n = 1; n <= part; n <<= 1)
		part -= n;

	double n2 = 2.0 / n;							// length of the part
	lo = ((double) part * n2 - 1.0);
	hi = ((double) (part + 1) * n2 - 1.0);
}

///
/// Get the position on the edge
/// @param part[in] ID of the position on the edge (see PIC)
/// @param x[out] position on the edge
void get_edge_part(int part, double &x) {
	if (part == 0)
		x = -1.0;
	else if (part == 1)
		x = 1.0;
	else {
		double lo, hi;
		get_interval_part(part - 2, lo, hi);
		x = (lo + hi) / 2.0;
	}
}

double get_edge_coef(int part) {
	// FIXME: handle endpoints (i.e. part == 0, part == 1) separately (?)
	double x;
	get_edge_part(part, x);
	return (x + 1.0) / 2.0;
}

Part transform_edge_part(int ori, Part part) {
	Part rp;
	rp.part = (ori == 0) ? part.part : opposite_part(part.part);
	return rp;
}


Part transform_face_part(int ori, Part part) {
	// refer to Pavel Solin's gray book, p. 169 (?)
	int flags[8][3] = {
		{ 1, 1, 1 }, { -1, 1, 1 }, { 1, -1, 1 }, { -1, -1, 1 }, { 1, 1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, -1 }
	};

	Part rp;
	if (flags[ori][2] == 1) {
		rp.horz = (flags[ori][0] > 0) ? part.horz : opposite_part(part.horz);
		rp.vert = (flags[ori][1] > 0) ? part.vert : opposite_part(part.vert);
	}
	else {
		// switch hpart and vpart
		rp.horz = (flags[ori][1] > 0) ? part.vert : opposite_part(part.vert);
		rp.vert = (flags[ori][0] > 0) ? part.horz : opposite_part(part.horz);
	}

	return rp;
}


// Shapeset /////

Shapeset::Shapeset() {
	mode = 0;
	ced_idx = -1;
	num_components = -1;

#ifdef PRELOADING
	fn_prods = NULL;
	dx_prods = NULL;
	dy_prods = NULL;
	dz_prods = NULL;
#endif
}

Shapeset::~Shapeset() {
#ifdef PRELOADING
	delete [] fn_prods;
	delete [] dx_prods;
	delete [] dy_prods;
	delete [] dz_prods;
#endif
	free_constrained_combinations();
}

int Shapeset::get_constrained_edge_index(int edge, int ori, order1_t order, Part part) {
	CEDKey cedkey(CED_KEY_TYPE_EDGE, edge, order, ori, part);
	int fn_idx;
	if (ced_id.lookup(cedkey, fn_idx))
		return -1 - fn_idx;
	else {
		// not existing yet => create an id for this situation
		ced_idx++;
		ced_key[ced_idx] = cedkey;
		ced_id.set(cedkey, ced_idx);
		return -1 - ced_idx;
	}
}

int Shapeset::get_constrained_edge_face_index(int edge, int ori, order2_t order, Part part, int dir) {
	CEDKey ck(CED_KEY_TYPE_EDGE_FACE, edge, order, ori, part, dir);
	int fn_idx;
	if (ced_id.lookup(ck, fn_idx))
		return -1 - fn_idx;
	else {
		// not existing yet => create an id for this situation
		ced_idx++;
		ced_key[ced_idx] = ck;
		ced_id.set(ck, ced_idx);
		return -1 - ced_idx;
	}
}

int Shapeset::get_constrained_face_index(int face, int ori, order2_t order, Part part) {
	CEDKey cedkey(CED_KEY_TYPE_FACE, face, order, ori, part);
	int fn_idx;
	if (ced_id.lookup(cedkey, fn_idx))
		return -1 - fn_idx;
	else {
		// not existing yet => create an id for this situation
		ced_idx++;
		ced_key[ced_idx] = cedkey;
		ced_id.set(cedkey, ced_idx);
		return -1 - ced_idx;
	}
}

void Shapeset::free_constrained_combinations() {
	for (Word_t i = ced_comb.first(); i != INVALID_IDX; i = ced_comb.next(i))
		delete ced_comb.get(i);
	ced_id.remove_all();
	ced_key.remove_all();
	ced_idx = -1;
}

CEDComb *Shapeset::get_ced_comb(const CEDKey &key) {
	CEDComb *comb;
	if (ced_comb.lookup(key, comb)) {
		// ok, already calculated combination
	}
	else {
		// combination does not exist yet => calculate it
		if (key.type == CED_KEY_TYPE_EDGE)           comb = calc_constrained_edge_combination(key.ori, key.order, key.part);
		else if (key.type == CED_KEY_TYPE_EDGE_FACE) comb = calc_constrained_edge_face_combination(key.ori, key.order, key.part, key.dir);
		else if (key.type == CED_KEY_TYPE_FACE)      comb = calc_constrained_face_combination(key.ori, key.order, key.part);
		else EXIT(ERR_FAILURE, "Unknown type of CED key.");

		ced_comb.set(key, comb);
	}

	return comb;
}

int *Shapeset::get_ced_indices(const CEDKey &key) {
	int *idx;
	if (key.type == CED_KEY_TYPE_EDGE) {
		order1_t order = key.order;
		idx = get_edge_indices(key.edge, key.ori, order);
	}
	else if (key.type == CED_KEY_TYPE_EDGE_FACE) {
		// HEX specific
		int dir = key.dir;
		const int *eori = RefHex::get_face_edge_orientation(key.ori);
		if (key.ori >= 4) dir = (key.dir == PART_ORI_VERT) ? PART_ORI_HORZ : PART_ORI_VERT; 			// turned face
		order2_t o = order2_t::from_int(key.order);
		idx = (dir == PART_ORI_HORZ) ? get_edge_indices(key.edge, eori[0], o.x) : get_edge_indices(key.edge, eori[1], o.y);
	}
	else if (key.type == CED_KEY_TYPE_FACE) {
		order2_t order = order2_t::from_int(key.order);
		idx = get_face_indices(key.face, key.ori, order);
	}
	else
		EXIT(ERR_FAILURE, "Unknown type of CED key.");

	return idx;
}

double Shapeset::get_constrained_value(int n, int index, double x, double y, double z, int component) {
	assert(ced_key.exists(-1 - index));
	CEDKey key = ced_key[-1 - index];

	CEDComb *comb = get_ced_comb(key);
	assert(comb != NULL);
	int *idx = get_ced_indices(key);
	assert(idx != NULL);

	double sum = 0.0;
	for (int i = 0; i < comb->n; i++)
		sum += comb->coef[i] * get_val(n, idx[i], x, y, z, component);

	return sum;
}

#ifdef PRELOADING

bool Shapeset::load_prods(const char *file_name, double *&mat) {
	FILE *file = fopen(file_name, "r");
	if (file != NULL) {
		fread(&num_fns, sizeof(num_fns), 1, file);
		int *idx_vec = new int [num_fns];
		MEM_CHECK(idx_vec);
		fread(idx_vec, sizeof(int), num_fns, file);

		if (fnidx2idx.count() == 0) {
			for (int i = 0; i < num_fns; i++)
				fnidx2idx[idx_vec[i]] = i;
		}

		mat = new double [num_fns * num_fns];
		MEM_CHECK(mat);
		fread(mat, sizeof(double), num_fns * num_fns, file);

		fclose(file);

		delete [] idx_vec;

		return true;
	}
	else
		return false;
}

bool Shapeset::preload_products() {
	// we do not need this for fichera
//	printf("Loading FN products.\n");
//	if (!load_prods("fn-fn", fn_prods)) return false;
	printf("Loading DX products.\n");
	if (!load_prods("dx-dx", dx_prods)) return false;
	printf("Loading DY products.\n");
	if (!load_prods("dy-dy", dy_prods)) return false;
	printf("Loading DZ products.\n");
	if (!load_prods("dz-dz", dz_prods)) return false;
	return true;
}

scalar Shapeset::get_product_val(int idx1, int idx2, double *vals) {
	assert(fnidx2idx.count() > 0 && vals != NULL);

	int idx[] = { idx1, idx2 };
	double ci[] = { 1.0, 1.0 };

	int nfn[2];
	int *fnidx[2];
	double *coef[2];
	for (int k = 0; k < 2; k++) {
		if (idx[k] < 0) {
			// ced
			assert(ced_key.exists(-1 - idx[k]));
			CEDKey key = ced_key[-1 - idx[k]];

			CEDComb *comb = get_ced_comb(key);

			fnidx[k] = get_ced_indices(key);
			nfn[k] = comb->n;
			coef[k] = comb->coef;
		}
		else {
			fnidx[k] = idx + k;
			nfn[k] = 1;
			coef[k] = ci + k;
		}
	}

	scalar res = 0.0;
	for (int i = 0; i < nfn[0]; i++) {
		scalar rj = 0.0;
		for (int j = 0; j < nfn[1]; j++)
			rj += coef[1][j] * vals[num_fns * fnidx2idx[fnidx[0][i]] + fnidx2idx[fnidx[1][j]]];
		res += coef[0][i] * rj;
	}

	return res;
}

#endif
