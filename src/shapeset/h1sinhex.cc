//
// h1sinhex.cc
//

#include "../config.h"
#include "common.h"
#include "h1sinhex.h"
#include <common/error.h>
#include "matrix.h"

#include "mesh.h"

#ifdef WITH_HEX

static const int base_coding = MAX_ELEMENT_ORDER + 1;		// order: 0..10

#define MAKE_HEX_IDX(h, f, v, ori) (((((v) * (base_coding) + (f)) * (base_coding) + (h)) << 3) | ori)

#define GET_HEX_IDX_1(order) ((order) % (base_coding))
#define GET_HEX_IDX_2(order) (((order) / (base_coding)) % (base_coding))
#define GET_HEX_IDX_3(order) ((order) / ((base_coding) * (base_coding)))

//#define MAKE_FACE_IDX(h, v) (((v) * (base_coding) + (h)) << 3)


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
/*
static scalar fn_2(double x)  { return sin((M_PI / 2) * (x + 1)); }
static scalar fn_3(double x)  { return sin(M_PI * (x + 1)); }
static scalar fn_4(double x)  { return sin(3 * (M_PI / 2) * (x + 1)); }
static scalar fn_5(double x)  { return sin(2 * M_PI * (x + 1)); }
static scalar fn_6(double x)  { return sin(5 * (M_PI / 2) * (x + 1)); }
static scalar fn_7(double x)  { return sin(3 * M_PI * (x + 1)); }
static scalar fn_8(double x)  { return sin(7 * (M_PI / 2) * (x + 1)); }
static scalar fn_9(double x)  { return sin(4 * M_PI * (x + 1)); }
static scalar fn_10(double x) { return sin(9 * (M_PI / 2) * (x + 1)); }
*/

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
/*
static scalar der_2(double x)  { return (M_PI / 2) * cos((M_PI / 2) * (x + 1)); }
static scalar der_3(double x)  { return M_PI * cos(M_PI * (x + 1)); }
static scalar der_4(double x)  { return 3 * (M_PI / 2) * cos(3 * (M_PI / 2) * (x + 1)); }
static scalar der_5(double x)  { return 2 * M_PI * cos(2 * M_PI * (x + 1)); }
static scalar der_6(double x)  { return 5 * (M_PI / 2) * cos(5 * (M_PI / 2) * (x + 1)); }
static scalar der_7(double x)  { return 3 * M_PI * cos(3 * M_PI * (x + 1)); }
static scalar der_8(double x)  { return 7 * (M_PI / 2) * cos(7 * (M_PI / 2) * (x + 1)); }
static scalar der_9(double x)  { return 4 * M_PI * cos(4 * M_PI * (x + 1)); }
static scalar der_10(double x) { return 9 * (M_PI / 2) * cos(9 * (M_PI / 2) * (x + 1)); }
*/

static shape_fn_1d_t der_tab_1d[] = { der_0, der_1, der_2, der_3, der_4, der_5, der_6, der_7, der_8, der_9, der_10 };

static int hex_vertex_indices[] = {
	MAKE_HEX_IDX(0, 0, 0, 0), MAKE_HEX_IDX(1, 0, 0, 0), MAKE_HEX_IDX(1, 1, 0, 0), MAKE_HEX_IDX(0, 1, 0, 0),
	MAKE_HEX_IDX(0, 0, 1, 0), MAKE_HEX_IDX(1, 0, 1, 0), MAKE_HEX_IDX(1, 1, 1, 0), MAKE_HEX_IDX(0, 1, 1, 0)
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


static void decompose(int index, int ori, int indices[3], int oris[3], bool swapori = true) {
	assert(ori >= 0);
	int permut[3];
	int num_01;

	// get order in each direction
	indices[0] = GET_HEX_IDX_1(index);
	indices[1] = GET_HEX_IDX_2(index);
	indices[2] = GET_HEX_IDX_3(index);

	find_permutation(indices, permut, num_01);

	oris[0] = oris[1] = oris[2] = 0;
	if (num_01 == 2) {
		// edge function
		assert(ori <= 1);
		oris[permut[2]] = ori;
	}
	else if (num_01 == 1) {
		// face function
		assert(ori <= 7);
		if (ori % 2 == 1) oris[permut[1]] = 1;
		if (ori % 4 >= 2) oris[permut[2]] = 1;
		if (ori >= 4) {
			std::swap(indices[permut[1]], indices[permut[2]]);
			if (swapori) std::swap(oris[permut[1]], oris[permut[2]]);
		}
	}
	else {
		// vertex or bubble function
		assert(ori == 0);
	}
}

// -- functions that calculate values of fn, dx, dy, dz on the fly -- //

static double calc_fn_value(int index, double x, double y, double z, int component) {
	int ori = GET_ORI_FROM_INDEX(index);
	int idx = GET_IDX_FROM_INDEX(index);

	int indices[3];
	int oris[3];

	decompose(idx, ori, indices, oris);

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	return fn_tab_1d[indices[0]](x) * fn_tab_1d[indices[1]](y) * fn_tab_1d[indices[2]](z);
}


static double calc_dx_value(int index, double x, double y, double z, int component) {
	int ori = GET_ORI_FROM_INDEX(index);
	int idx = GET_IDX_FROM_INDEX(index);

	int indices[3];
	int oris[3];

	decompose(idx, ori, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	double dx = der_tab_1d[indices[0]](x) * fn_tab_1d[indices[1]](y) * fn_tab_1d[indices[2]](z);
	if (oris[0] == 1) dx = -dx;

	return dx;
}


static double calc_dy_value(int index, double x, double y, double z, int component) {
	int ori = GET_ORI_FROM_INDEX(index);
	int idx = GET_IDX_FROM_INDEX(index);

	int indices[3];
	int oris[3];

	decompose(idx, ori, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	double dy = fn_tab_1d[indices[0]](x) * der_tab_1d[indices[1]](y) * fn_tab_1d[indices[2]](z);
	if (oris[1] == 1) dy = -dy;

	return dy;
}


static double calc_dz_value(int index, double x, double y, double z, int component) {
	int ori = GET_ORI_FROM_INDEX(index);
	int idx = GET_IDX_FROM_INDEX(index);

	int indices[3];
	int oris[3];

	decompose(idx, ori, indices, oris);

	for (int i = 0; i < 3; i++)
		assert((oris[i] == 0) || (indices[i] >= 2));

	if (oris[0] == 1) x = -x;
	if (oris[1] == 1) y = -y;
	if (oris[2] == 1) z = -z;

	double dz = fn_tab_1d[indices[0]](x) * fn_tab_1d[indices[1]](y) * der_tab_1d[indices[2]](z);
	if (oris[2] == 1) dz = -dz;

	return dz;
}

#endif

// -- -- //

H1ShapesetSinHex::H1ShapesetSinHex() {
#ifdef WITH_HEX
	mode = MODE_HEXAHEDRON;

	num_components = 1;

	max_edge_order = MAX_ELEMENT_ORDER + 1;
	max_face_order = MAKE_QUAD_ORDER(MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1);
	max_order = MAKE_HEX_ORDER(MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1);

	max_index = MAKE_HEX_IDX(MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1, MAX_ELEMENT_ORDER + 1, 0);

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

	// index to order mapping
	index_to_order = new int[max_order];
	MEM_CHECK(index_to_order);
	for (int i = 0; i <= MAX_ELEMENT_ORDER; i++)
		for (int j = 0; j <= MAX_ELEMENT_ORDER; j++)
			for (int k = 0; k <= MAX_ELEMENT_ORDER; k++)
				index_to_order[MAKE_HEX_IDX(i, j, k, 0) >> 3] = MAKE_HEX_ORDER(index_order[i], index_order[j], index_order[k]);
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

H1ShapesetSinHex::~H1ShapesetSinHex() {
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

	delete [] index_to_order;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

void H1ShapesetSinHex::compute_edge_indices(int edge, int ori, int order) {
#ifdef WITH_HEX
	int *indices = new int[order - 1];
	MEM_CHECK(indices);

	int idx = 0;
	switch (edge) {
		case  0: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(i, 0, 0, ori); break;
		case  1: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(1, i, 0, ori); break;
		case  2: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(i, 1, 0, ori); break;
		case  3: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(0, i, 0, ori); break;
		case  4: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(0, 0, i, ori); break;
		case  5: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(1, 0, i, ori); break;
		case  6: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(1, 1, i, ori); break;
		case  7: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(0, 1, i, ori); break;
		case  8: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(i, 0, 1, ori); break;
		case  9: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(1, i, 1, ori); break;
		case 10: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(i, 1, 1, ori); break;
		case 11: for (int i = 2; i <= order; i++) indices[idx++] = MAKE_HEX_IDX(0, i, 1, ori); break;
		default: EXIT(ERR_FAILURE, "Invalid edge number %d. Can be 0 - 11.", edge); break;
	}

	edge_indices[edge][ori][order] = indices;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

void H1ShapesetSinHex::compute_face_indices(int face, int ori, int order) {
#ifdef WITH_HEX
	int *indices = new int[(GET_QUAD_ORDER_1(order) - 1) * (GET_QUAD_ORDER_2(order) - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	switch (face) {
		case 0:
			for(int i = 2; i <= GET_QUAD_ORDER_1(order); i++)
				for(int j = 2; j <= GET_QUAD_ORDER_2(order); j++)
					indices[idx++] = MAKE_HEX_IDX(0, i, j, ori);
			break;

		case 1:
			for(int i = 2; i <= GET_QUAD_ORDER_1(order); i++)
				for(int j = 2; j <= GET_QUAD_ORDER_2(order); j++)
					indices[idx++] = MAKE_HEX_IDX(1, i, j, ori);
			break;

		case 2:
			for(int i = 2; i <= GET_QUAD_ORDER_1(order); i++)
				for(int j = 2; j <= GET_QUAD_ORDER_2(order); j++)
					indices[idx++] = MAKE_HEX_IDX(i, 0, j, ori);
			break;

		case 3:
			for(int i = 2; i <= GET_QUAD_ORDER_1(order); i++)
				for(int j = 2; j <= GET_QUAD_ORDER_2(order); j++)
					indices[idx++] = MAKE_HEX_IDX(i, 1, j, ori);
			break;

		case 4:
			for(int i = 2; i <= GET_QUAD_ORDER_1(order); i++)
				for(int j = 2; j <= GET_QUAD_ORDER_2(order); j++)
					indices[idx++] = MAKE_HEX_IDX(i, j, 0, ori);
			break;

		case 5:
			for(int i = 2; i <= GET_QUAD_ORDER_1(order); i++)
				for(int j = 2; j <= GET_QUAD_ORDER_2(order); j++)
					indices[idx++] = MAKE_HEX_IDX(i, j, 1, ori);
			break;

		default:
			EXIT(ERR_FAILURE, "Invalid face number %d. Can be 0 - 5.", face);
			break;
	}

	face_indices[face][ori][order] = indices;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

void H1ShapesetSinHex::compute_bubble_indices(int order) {
#ifdef WITH_HEX
	int order1 = GET_HEX_ORDER_1(order);
	int order2 = GET_HEX_ORDER_2(order);
	int order3 = GET_HEX_ORDER_3(order);
	int *indices = new int[(order1 - 1) * (order2 - 1) * (order3 - 1)];
	MEM_CHECK(indices);

	int idx = 0;
	for (int i = 2; i <= order1; i++)
		for(int j = 2; j <= order2; j++)
			for(int k = 2; k <= order3; k++)
				indices[idx++] = MAKE_HEX_IDX(i, j, k, 0);

	bubble_indices[order] = indices;
#else
	EXIT(ERR_HEX_NOT_COMPILED);
#endif
}

/// --- CED specific stuff ---

CEDComb *H1ShapesetSinHex::calc_constrained_edge_combination(int ori, int order, Part part) {
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}

CEDComb *H1ShapesetSinHex::calc_constrained_edge_face_combination(int ori, int order, Part part) {
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}

CEDComb *H1ShapesetSinHex::calc_constrained_face_combination(int ori, int order, Part part) {
	EXIT(ERR_NOT_IMPLEMENTED);
	return NULL;
}

