//
// h1lobbatohex.cc
//

#include "../config.h"
#include "common.h"
#include "lobatto.h"
#include "h1lobattohex.h"
#include <common/error.h>
#include "matrix.h"

#include "mesh.h"

#ifdef WITH_HEX

struct hex_index_t {
	unsigned ori:3;
	unsigned x:8;
	unsigned y:8;
	unsigned z:8;

	hex_index_t(int idx) {
		this->x = (idx >> 16) & 0x0F;
		this->y = (idx >>  8) & 0x0F;
		this->z = (idx >>  0) & 0x0F;
		this->ori = (idx >> 24) & 0x07;
	}

	hex_index_t(int x, int y, int z, int ori = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->ori = ori;
	}

	operator int() { return (((((ori << 8) | x) << 8) | y) << 8) | z; }
};


static int lobatto_hex_vertex_indices[] = {
	hex_index_t(0, 0, 0), hex_index_t(1, 0, 0), hex_index_t(1, 1, 0), hex_index_t(0, 1, 0),
	hex_index_t(0, 0, 1), hex_index_t(1, 0, 1), hex_index_t(1, 1, 1), hex_index_t(0, 1, 1)
};

static int index_order[] = { 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

// -- helpers -- //

static void find_permutation(int *indices, int *permut, int &num_01) {
	for (int i = 0; i < 3; i++)
		permut[i] = i;
	num_01 = 0;
	if (indices[0] < 2) num_01++;
	if (indices[1] < 2) {
		num_01++;
		if (num_01 == 1) swapint(permut[0], permut[1]);
	}
	if (indices[2] < 2) {
		num_01++;
		if (num_01 == 1) {
			swapint(permut[1], permut[2]);
			swapint(permut[0], permut[1]);
		}
		if (num_01 == 2) swapint(permut[1], permut[2]);
	}
}


static void decompose(hex_index_t index, int indices[3], int ori[3], bool swapori = true) {
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
			swapint(indices[permut[1]], indices[permut[2]]);
			if (swapori) swapint(ori[permut[1]], ori[permut[2]]);
		}
	}
	else {
		// vertex or bubble function
		assert(index.ori == 0);
	}
}

// -- functions that calculate values of fn, dx, dy, dz on the fly -- //

static double calc_fn_value(int index, double x, double y, double z, int component) {
	hex_index_t idx(index);
	int indices[3];
	int oris[3];

	decompose(idx, indices, oris);

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	return lobatto_fn_tab_1d[indices[0]](x) * lobatto_fn_tab_1d[indices[1]](y) * lobatto_fn_tab_1d[indices[2]](z);
}


static double calc_dx_value(int index, double x, double y, double z, int component) {
	hex_index_t idx(index);
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
	hex_index_t idx(index);
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
	hex_index_t idx(index);
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

H1ShapesetLobattoHex::H1ShapesetLobattoHex() {
#ifdef WITH_HEX
	mode = MODE_HEXAHEDRON;

	num_components = 1;

//	max_edge_order = MAX_ELEMENT_ORDER + 1;
//	max_face_order = MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1);
//	max_order = MAKE_HEX_ORDER(MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1);
//	max_index = MAKE_HEX_IDX(MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1, 0);

	// fn, dx, dy, dz will be calculated on-the-fly
	shape_table_deleg[FN]  = calc_fn_value;
	shape_table_deleg[DX]  = calc_dx_value;
	shape_table_deleg[DY]  = calc_dy_value;
	shape_table_deleg[DZ]  = calc_dz_value;
	shape_table_deleg[DXY] = NULL;
	shape_table_deleg[DXZ] = NULL;
	shape_table_deleg[DYZ] = NULL;

	// vertices
	vertex_indices = lobatto_hex_vertex_indices;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

H1ShapesetLobattoHex::~H1ShapesetLobattoHex() {
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
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

order3_t H1ShapesetLobattoHex::get_order(int index) const {
	if (index >= 0) {
		hex_index_t idx(index);
		order3_t ord(index_order[idx.x], index_order[idx.y], index_order[idx.z]);
		if (idx.ori >= 4) ord = turn_hex_face_order(ord);		// face function is turned due to orientation
		return ord;
	}
	else {
		assert(ced_key.exists(-1 - index));
		CEDKey key = ced_key[-1 - index];

		int order;
		if (key.type == CED_KEY_TYPE_EDGE || key.type == CED_KEY_TYPE_EDGE_FACE) {
//			if (key.type == CED_KEY_TYPE_EDGE_FACE) order = (key.ori < 4) ? GET_QUAD_ORDER_1(key.order) : GET_QUAD_ORDER_2(key.order);
//			else order = key.order;

//			switch (key.edge) {
//				case  0: order = MAKE_HEX_ORDER(order, 1, 1); break;
//				case  1: order = MAKE_HEX_ORDER(1, order, 1); break;
//				case  2: order = MAKE_HEX_ORDER(order, 1, 1); break;
//				case  3: order = MAKE_HEX_ORDER(1, order, 1); break;
//				case  4: order = MAKE_HEX_ORDER(1, 1, order); break;
//				case  5: order = MAKE_HEX_ORDER(1, 1, order); break;
//				case  6: order = MAKE_HEX_ORDER(1, 1, order); break;
//				case  7: order = MAKE_HEX_ORDER(1, 1, order); break;
//				case  8: order = MAKE_HEX_ORDER(order, 1, 1); break;
//				case  9: order = MAKE_HEX_ORDER(1, order, 1); break;
//				case 10: order = MAKE_HEX_ORDER(order, 1, 1); break;
//				case 11: order = MAKE_HEX_ORDER(1, order, 1); break;
//				default: EXIT(ERR_FAILURE, "Invalid edge number %d. Can be 0 - 11.", key.edge); break;
//			}
		}
		else if (key.type == CED_KEY_TYPE_FACE) {
//			int o[] = { GET_QUAD_ORDER_1(key.order), GET_QUAD_ORDER_2(key.order) };
//
//			switch (key.face) {
//				case 0: order = MAKE_HEX_ORDER(1, o[0], o[1]); break;
//				case 1: order = MAKE_HEX_ORDER(1, o[0], o[1]); break;
//				case 2: order = MAKE_HEX_ORDER(o[0], 1, o[1]); break;
//				case 3: order = MAKE_HEX_ORDER(o[0], 1, o[1]); break;
//				case 4: order = MAKE_HEX_ORDER(o[0], o[1], 1); break;
//				case 5: order = MAKE_HEX_ORDER(o[0], o[1], 1); break;
//				default: EXIT(ERR_FAILURE, "Invalid face number %d. Can be 0 - 5.", key.face); break;
//			}
		}
		if (key.ori >= 4) return turn_hex_face_order(order);
		else return order;
	}
}

void H1ShapesetLobattoHex::compute_edge_indices(int edge, int ori, order1_t order) {
#ifdef WITH_HEX
	int *indices = new int[order - 1];
	MEM_CHECK(indices);

	int idx = 0;
	switch (edge) {
		case  0: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(i, 0, 0, ori); break;
		case  1: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(1, i, 0, ori); break;
		case  2: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(i, 1, 0, ori); break;
		case  3: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(0, i, 0, ori); break;
		case  4: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(0, 0, i, ori); break;
		case  5: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(1, 0, i, ori); break;
		case  6: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(1, 1, i, ori); break;
		case  7: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(0, 1, i, ori); break;
		case  8: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(i, 0, 1, ori); break;
		case  9: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(1, i, 1, ori); break;
		case 10: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(i, 1, 1, ori); break;
		case 11: for (int i = 2; i <= order; i++) indices[idx++] = hex_index_t(0, i, 1, ori); break;
		default: EXIT(ERR_FAILURE, "Invalid edge number %d. Can be 0 - 11.", edge); break;
	}

	edge_indices[edge][ori][order] = indices;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

void H1ShapesetLobattoHex::compute_face_indices(int face, int ori, order2_t order) {
#ifdef WITH_HEX
	int horder = order.x, vorder = order.y;
	int *indices = new int[(horder - 1) * (vorder - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	switch (face) {
		case 0:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = hex_index_t(0, i, j, ori);
			break;

		case 1:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = hex_index_t(1, i, j, ori);
			break;

		case 2:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = hex_index_t(i, 0, j, ori);
			break;

		case 3:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = hex_index_t(i, 1, j, ori);
			break;

		case 4:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = hex_index_t(i, j, 0, ori);
			break;

		case 5:
			for(int i = 2; i <= horder; i++)
				for(int j = 2; j <= vorder; j++)
					indices[idx++] = hex_index_t(i, j, 1, ori);
			break;

		default:
			EXIT(ERR_FAILURE, "Invalid face number %d. Can be 0 - 5.", face);
			break;
	}

	face_indices[face][ori][order.get_idx()] = indices;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

void H1ShapesetLobattoHex::compute_bubble_indices(order3_t order) {
#ifdef WITH_HEX
	int *indices = new int[(order.x - 1) * (order.y - 1) * (order.z - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	for (int i = 2; i <= order.x; i++)
		for(int j = 2; j <= order.y; j++)
			for(int k = 2; k <= order.z; k++)
				indices[idx++] = hex_index_t(i, j, k, 0);

	bubble_indices[order.get_idx()] = indices;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

/// --- CED specific stuff ---

static Part transform_edge_part(int ori, Part part) {
	Part rp;
	rp.part = (ori == 0) ? part.part : opposite_part(part.part);
	return rp;
}

//
// constraints are calculated on egde 0
//
CEDComb *H1ShapesetLobattoHex::calc_constrained_edge_combination(int ori, order1_t order, Part part) {
#ifdef WITH_HEX
	Part rp = transform_edge_part(ori, part);

	// determine the interval of the edge
	double hi, lo;
	get_interval_part(rp.part, lo, hi);

	int n = get_num_edge_fns(order);							// total number of functions on the edge
	int *fn_idx = get_edge_indices(0, 0, order);				// indices of all functions on the edge

	double f_lo = get_value(FN, fn_idx[n - 1], lo, -1.0, -1.0, 0);		// fn. values at endpoints of the part
	double f_hi = get_value(FN, fn_idx[n - 1], hi, -1.0, -1.0, 0);

	double **a = new_matrix<double>(n, n);
	MEM_CHECK(a);
	double *b = new double[n];
	MEM_CHECK(b);
	for (int i = 0; i < n; i++) {
		// chebyshev point
		double p = cos((i + 1) * M_PI / order);
		double r = (p + 1.0) * 0.5;
		double s = 1.0 - r;

		// matrix row
		for (int j = 0; j < n; j++)
			a[i][j] = get_value(FN, fn_idx[j], p, -1.0, -1.0, 0);

		// rhs
		b[i] = get_value(FN, fn_idx[n - 1], lo*s + hi*r, -1.0, -1.0, 0) - f_lo*s - f_hi*r;
	}

	// solve the system
	double d;
	int *iperm = new int[n];
	MEM_CHECK(iperm);
	ludcmp(a, n, iperm, &d);
	lubksb(a, n, iperm, b);

	// cleanup
	delete [] iperm;
	delete [] a;

	return new CEDComb(n, b);
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

static Part transform_face_part(int ori, Part part) {
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

//
// constraints are calculated on face 5
//
CEDComb *H1ShapesetLobattoHex::calc_constrained_edge_face_combination(int ori, order2_t order, Part part, int dir) {
#ifdef WITH_HEX
	Part rp = transform_face_part(ori, part);

	if (ori >= 4) dir = (dir == PART_ORI_VERT) ? PART_ORI_HORZ : PART_ORI_VERT; 			// turned face

	// determine the interval of the face
	double hi, lo;
	double x0;
	int epart;

	if (dir == PART_ORI_VERT) {
		get_interval_part(rp.vert, lo, hi);
		epart = face_to_edge_part(rp.horz);
		get_edge_part(epart, x0);

		int horder = 0;//GET_QUAD_ORDER_1(order);
		int vorder = 0;//GET_QUAD_ORDER_2(order);

		int n = get_num_edge_fns(vorder);										// total number of functions on the edge
		int *edge_fn_idx[] = {
			get_edge_indices(0, 0, horder),										// indices of all functions on the edge
			get_edge_indices(0, 0, vorder)										// indices of all functions on the edge
		};

		double f_lo = get_value(FN, edge_fn_idx[1][n - 1], lo, -1.0, -1.0, 0);		// fn. values at endpoints of the part
		double f_hi = get_value(FN, edge_fn_idx[1][n - 1], hi, -1.0, -1.0, 0);

		double **a = new_matrix<double>(n, n);
		MEM_CHECK(a);
		double *b = new double[n];
		MEM_CHECK(b);
		for (int i = 0; i < n; i++) {
			// chebyshev point
			double p = cos((i + 1) * M_PI / vorder);
			double r = (p + 1.0) * 0.5;
			double s = 1.0 - r;

			// matrix row
			for (int j = 0; j < n; j++)
				a[i][j] = get_value(FN, edge_fn_idx[1][j], p, -1.0, -1.0, 0);

			// rhs
			b[i] = get_value(FN, edge_fn_idx[1][n - 1], lo*s + hi*r, -1.0, -1.0, 0) - f_lo*s - f_hi*r;	// depends on the ref. domain
		}

		// solve the system
		double d;
		int *iperm = new int[n];
		MEM_CHECK(iperm);
		ludcmp(a, n, iperm, &d);
		lubksb(a, n, iperm, b);

		int m = get_num_edge_fns(horder);											// total number of functions on the edge
		double c = get_value(FN, edge_fn_idx[0][m - 1], x0, -1.0, -1.0, 0);
		for (int i = 0; i < n; i++)
			b[i] *= c;

		// cleanup
		delete [] iperm;
		delete [] a;

		return new CEDComb(n, b);
	}
	else {
		get_interval_part(rp.horz, lo, hi);
		epart = face_to_edge_part(rp.vert);
		get_edge_part(epart, x0);

		int horder = 0;//GET_QUAD_ORDER_1(order);
		int vorder = 0;//GET_QUAD_ORDER_2(order);

		int n = get_num_edge_fns(horder);										// total number of functions on the edge
		int *edge_fn_idx[] = {
			get_edge_indices(0, 0, horder),										// indices of all functions on the edge
			get_edge_indices(0, 0, vorder)										// indices of all functions on the edge
		};

		double f_lo = get_value(FN, edge_fn_idx[0][n - 1], lo, -1.0, -1.0, 0);		// fn. values at endpoints of the part
		double f_hi = get_value(FN, edge_fn_idx[0][n - 1], hi, -1.0, -1.0, 0);

		double **a = new_matrix<double>(n, n);
		MEM_CHECK(a);
		double *b = new double[n];
		MEM_CHECK(b);
		for (int i = 0; i < n; i++) {
			// chebyshev point
			double p = cos((i+1) * M_PI / horder);
			double r = (p + 1.0) * 0.5;
			double s = 1.0 - r;

			// matrix row
			for (int j = 0; j < n; j++)
				a[i][j] = get_value(FN, edge_fn_idx[0][j], p, -1.0, -1.0, 0);

			// rhs
			b[i] = get_value(FN, edge_fn_idx[0][n - 1], lo*s + hi*r, -1.0, -1.0, 0) - f_lo*s - f_hi*r;
		}

		// solve the system
		double d;
		int *iperm = new int[n];
		MEM_CHECK(iperm);
		ludcmp(a, n, iperm, &d);
		lubksb(a, n, iperm, b);

		int m = get_num_edge_fns(vorder);											// total number of functions on the edge
		double c = get_value(FN, edge_fn_idx[1][m - 1], x0, -1.0, -1.0, 0);
		for (int i = 0; i < n; i++)
			b[i] *= c;

		// cleanup
		delete [] iperm;
		delete [] a;

		return new CEDComb(n, b);
	}
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

//
// constraints are calculated on face 5
//
//
//            edge2
//  v_hi +-----------+
//       |           |
// edge3 |           | edge1
//       |           |
//  v_lo +-----------+
//     h_lo  edge0  h_hi
//
CEDComb *H1ShapesetLobattoHex::calc_constrained_face_combination(int ori, order2_t order, Part part) {
#ifdef WITH_HEX
	order2_t old_order = order;

	int n = get_num_face_fns(order);										// total number of functions on the face
	int *fn_idx = get_face_indices(5, 0, order);							// indices of all functions on the face

	Part rp = transform_face_part(ori, part);

	double h_hi, h_lo, v_hi, v_lo;
	get_interval_part(rp.horz, h_lo, h_hi);				// determine the horizontal interval of the face
	get_interval_part(rp.vert, v_lo, v_hi);				// determine the vertical interval of the face

	int horder, vorder;
//	horder = GET_QUAD_ORDER_1(order);
//	vorder = GET_QUAD_ORDER_2(order);

	get_interval_part(rp.horz, h_lo, h_hi);
	get_interval_part(rp.vert, v_lo, v_hi);

	// fn. values at vertices of the face part
	double f_lo_lo = get_value(FN, fn_idx[n - 1], h_lo, v_lo, 1.0, 0);
	double f_lo_hi = get_value(FN, fn_idx[n - 1], h_lo, v_hi, 1.0, 0);
	double f_hi_lo = get_value(FN, fn_idx[n - 1], h_hi, v_lo, 1.0, 0);
	double f_hi_hi = get_value(FN, fn_idx[n - 1], h_hi, v_hi, 1.0, 0);

	// edge parts
	int *edge_fn_idx[4];
	edge_fn_idx[0] = get_edge_indices( 9, 0, vorder);
	edge_fn_idx[1] = get_edge_indices(10, 0, horder);
	edge_fn_idx[2] = get_edge_indices(11, 0, vorder);
	edge_fn_idx[3] = get_edge_indices( 8, 0, horder);

	double f_edge[4];
	f_edge[0] = get_fn_value(edge_fn_idx[0][vorder - 2],  1.0, v_lo, 1.0, 0);
	f_edge[1] = get_fn_value(edge_fn_idx[1][horder - 2], h_hi,  1.0, 1.0, 0);
	f_edge[2] = get_fn_value(edge_fn_idx[2][vorder - 2], -1.0, v_hi, 1.0, 0);
	f_edge[3] = get_fn_value(edge_fn_idx[3][horder - 2], h_lo, -1.0, 1.0, 0);

	Part hpart;
	Part vpart;
	hpart.part = rp.horz;
	vpart.part = rp.vert;
	int ced_edge_idx[4];
	ced_edge_idx[0] = get_constrained_edge_index( 8, 0, horder, hpart);
	ced_edge_idx[1] = get_constrained_edge_index( 9, 0, vorder, vpart);
	ced_edge_idx[2] = get_constrained_edge_index(10, 0, horder, hpart);
	ced_edge_idx[3] = get_constrained_edge_index(11, 0, vorder, vpart);

	double **a = new_matrix<double>(n, n);
	MEM_CHECK(a);
	double *b = new double[n];
	MEM_CHECK(b);

	//
	for (int row = 0; row < n; row++) {
		int face_order = 0;//get_hex_face_order(5, get_order(fn_idx[row]));
		int i = 0;//GET_QUAD_ORDER_1(face_order);
		int j = 0;//GET_QUAD_ORDER_2(face_order);

		double hp = cos((i - 1) * M_PI / horder);
		double hr = (hp + 1.0) * 0.5;
		double hs = 1.0 - hr;

		double vp = cos((j - 1) * M_PI / vorder);
		double vr = (vp + 1.0) * 0.5;
		double vs = 1.0 - vr;

		for (int k = 0; k < n; k++)
			a[row][k] = get_value(FN, fn_idx[k], hp, vp, 1.0, 0);

		// rhs
		b[row] = get_value(FN, fn_idx[n - 1], h_lo * hs + h_hi * hr, v_lo * vs + v_hi * vr, 1.0, 0)
			// subtract linear part
			- f_lo_lo * hs * vs
			- f_lo_hi * hs * vr
			- f_hi_lo * hr * vs
			- f_hi_hi * hr * vr
			// subtract residual of edge functions
			- f_edge[0] * get_constrained_value(FN, ced_edge_idx[0], hp, -1.0, 1.0, 0) * vs
			- f_edge[1] * get_constrained_value(FN, ced_edge_idx[1],  1.0, vp, 1.0, 0) * hr
			- f_edge[2] * get_constrained_value(FN, ced_edge_idx[2], hp,  1.0, 1.0, 0) * vr
			- f_edge[3] * get_constrained_value(FN, ced_edge_idx[3], -1.0, vp, 1.0, 0) * hs;
	}

	// solve the system
	double d;
	int *iperm = new int[n];
	MEM_CHECK(iperm);
	ludcmp(a, n, iperm, &d);
	lubksb(a, n, iperm, b);

	// cleanup
	delete [] iperm;
	delete [] a;

	return new CEDComb(n, b);
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}
