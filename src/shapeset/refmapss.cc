//
// refmapshapeset.cc
//
//

#include "../config.h"
#include "../shapeset.h"
#include "common.h"
#include "refmapss.h"
#include <common/error.h>
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

RefMapShapesetTetra::RefMapShapesetTetra() {
#ifdef WITH_TETRA
	mode = MODE_TETRAHEDRON;

	max_order = 10;			// MAX_ELEMENT_ORDER ?
	num_components = 1;

	max_index = (sizeof(refmap_tetra_fn) / sizeof(shape_fn_t)) - 1;

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
#ifdef WITH_TETRA
#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}

//// RefMapShapesetHex ////////////////////////////////////////////////////////////////////////////
#ifdef WITH_HEX

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
static int refmap_hex_index_to_order[] = { 1, 1, 1, 1, 1, 1, 1, 1 };

static shape_fn_t *refmap_hex_fn_table[] = { refmap_hex_fn };
static shape_fn_t *refmap_hex_dx_table[] = { refmap_hex_dx };
static shape_fn_t *refmap_hex_dy_table[] = { refmap_hex_dy };
static shape_fn_t *refmap_hex_dz_table[] = { refmap_hex_dz };

#endif


RefMapShapesetHex::RefMapShapesetHex() {
#ifdef WITH_HEX
	mode = MODE_HEXAHEDRON;

	max_order = 10;			// MAX_ELEMENT_ORDER ?
	num_components = 1;

	max_index = (sizeof(refmap_hex_fn) / sizeof(shape_fn_t)) - 1;

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
#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}

RefMapShapesetHex::~RefMapShapesetHex() {
#ifdef WITH_HEX
#else
	EXIT(ERR_TETRA_NOT_COMPILED);
#endif
}


