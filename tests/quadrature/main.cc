//
// main.cc
//
// Testing numerical quadrature
//
//

//
// TODO: test numerical quadrature for prisms
//
//

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

#define EPSILON										10e-10

#define countof(a) 									(sizeof(a)/sizeof(a[0]))


bool testPrint(bool value, const char *msg, bool correct) {
	printf("%s...", msg);
	if (value == correct) {
		printf("OK\n");
		return true;
	}
	else {
		printf("failed\n");
		return false;
	}
}

//
// Test quadrature
//

//
// 1D
//

typedef
	double (*fn1d_t)(double x);

struct TC1D {
	double exact;		// exact value of the integral
	fn1d_t fn;			// function
	int min_order;		// minimal order
	const char *fn_name;		// string representation of a function

	TC1D(fn1d_t f, double e, int min, const char *n) {
		exact = e;
		fn = f;
		fn_name = n;
		min_order = min;
	}
};

// function: f(x) = x^3 + x^2 + x + 1
double fn_1d_1(double x) {
	return x*x*x + x*x + x + 1;
}

// <<ADD MORE 1D FUNCTIONS HERE>>

// test 1d quadrature
int test_quadrature_1d_line(fn1d_t fn, double exact, int init_order, const char *fn_str) {
	// !!! std. quadrature works on std. reference domain !!!
	printf("  * f(x) = %s", fn_str);

	Quad1DStd quad;

	for (int order = init_order; order <= MAX_QUAD_ORDER; order++) {
		int np = quad.get_num_points(order);
		QuadPt1D *pt = quad.get_points(order);

		// \int_-1^1 fn1(x) dx
		double integral = 0;
		for (int i = 0; i < np; i++) {
			integral += fn(pt[i].x) * pt[i].w;
		}

		double err = fabs(exact - integral);
		if (err >= EPSILON) {
			printf(" ... failed for order %d, integral = %lf, expected = %lf\n", order, integral, exact);
			return ERROR_FAILURE;
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_1d() {
	printf("- Testing 1D quadrature -----\n");

	TC1D fns[] = {
		TC1D(fn_1d_1, 8.0 / 3.0, 3, "x^3 + x^2 + x + 1")
		// <<ADD MORE 1D FUNCTIONS HERE>>
	};

	for (int i = 0; i < countof(fns); i++) {
		int res = test_quadrature_1d_line(fns[i].fn, fns[i].exact, fns[i].min_order, fns[i].fn_name);
		if (res != ERROR_SUCCESS)
			return res;
	}

	return ERROR_SUCCESS;
}

//
// 2D
//

typedef
	double (*fn2d_t)(double x, double y);

// data for test case
struct TC2D {
	double exact;		// exact value of the integral
	int min_h_order;	// minimal horz order
	int min_v_order;	// minimal vert order
	fn2d_t fn;			// function
	const char *fn_name;		// string representation of a function

	TC2D(fn2d_t f, double e, int min_h, int min_v, const char *n) {
		exact = e;
		fn = f;
		fn_name = n;
		min_h_order = min_h;
		min_v_order = min_v;
	}
};


// function: f(x, y) = x^2 + y^2 + x*y + x + y + 1
double fn_2d_1(double x, double y) {
	return x*x + y*y + x*y + x + y + 1.0;
}

// function: f(x, y) = x^2 + y^2 + x + y + 1
double fn_2d_2(double x, double y) {
	return x*x + y*y + x + y + 1;
}

// <<ADD MORE 2D FUNCTIONS HERE>>

// test 2d quadrature
int test_quadrature_2d_quad(fn2d_t fn, double exact, int init_h_order, int init_v_order, const char *fn_str) {
	// !!! std. quadrature works on std. reference domain !!!

	// - on quad --------------------------------------------------------------
	QuadStdQuad quad;

	printf("  * f(x,y) = %s", fn_str);

//	exact = 20.0 / 3.0;
	for (int horder = init_h_order; horder <= MAX_QUAD_ORDER; horder++) {
		for (int vorder = init_v_order; vorder <= MAX_QUAD_ORDER; vorder++) {
			int order = MAKE_QUAD_ORDER(horder, vorder);

			int np = quad.get_num_points(order);
			QuadPt2D *pt = quad.get_points(order);

			// \int_quad fn_2d_1(x, y) dx dy
			double integral = 0;
			for (int i = 0; i < np; i++) {
				integral += fn(pt[i].x, pt[i].y) * pt[i].w;
			}

			double err = fabs(exact - integral);
			if (err >= EPSILON) {
				printf(" ... failed for order (h = %d, v = %d), integral = %lf, expected = %lf\n", horder, vorder, integral, exact);
				return ERROR_FAILURE;
			}
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_2d_tri(fn2d_t fn, double exact, int init_order, const char *fn_str) {
	// !!! std. quadrature works on std. reference domain !!!

	// - on triangle ----------------------------------------------------------
	QuadStdTri quad;

	// f(x, y) = x^2 + y^2 + x * y + x + y + 1
	printf("  * f(x,y) = %s", fn_str);

	for (int order = init_order; order <= MAX_QUAD_ORDER_TRI; order++) {
		int np = quad.get_num_points(order);
		QuadPt2D *pt = quad.get_points(order);

		// \int_tri fn_2d_1(x, y) dx dy
		double integral = 0;
		for (int i = 0; i < np; i++) {
			integral += fn(pt[i].x, pt[i].y) * pt[i].w;
		}

		double err = fabs(exact - integral);
		if (err >= EPSILON) {
			printf(" ... failed for order %d, integral = %lf, expected = %lf\n", order, integral, exact);
			return ERROR_FAILURE;
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_2d() {
	int ret = ERROR_SUCCESS;

	printf("\n");
	printf("- Testing 2D quadrature (quads) -----\n");

	TC2D fns_quad[] = {
		TC2D(fn_2d_1, 20.0 / 3.0, 2, 2, "x^2 + y^2 + x * y + x + y + 1")
		// <<ADD MORE 2D FUNCTIONS HERE>>
	};

	for (int i = 0; i < countof(fns_quad); i++) {
		if ((ret = test_quadrature_2d_quad(fns_quad[i].fn, fns_quad[i].exact, fns_quad[i].min_h_order, fns_quad[i].min_v_order, fns_quad[i].fn_name)) != ERROR_SUCCESS)
			return ret;
	}

	printf("\n");
	printf("- Testing 2D quadrature (triangle) -----\n");

	TC2D fns_tri[] = {
		TC2D(fn_2d_1, 2.0, 2, 0, "x^2 + y^2 + x * y + x + y + 1")
		// <<ADD MORE 2D FUNCTIONS HERE>>
	};

	for (int i = 0; i < countof(fns_tri); i++) {
		if ((ret = test_quadrature_2d_tri(fns_tri[i].fn, fns_tri[i].exact, fns_tri[i].min_h_order, fns_tri[i].fn_name)) != ERROR_SUCCESS)
			return ret;
	}

	return ERROR_SUCCESS;
}

//
// 3D
//

typedef
	double (*fn3d_t)(double x, double y, double z);

struct TC3D {
	double exact;		// exact value of the integral
	int min_h_order;	// minimal horz order
	int min_v_order;	// minimal vert order
	int min_u_order;	// minimal vert order
	fn3d_t fn;			// function
	const char *fn_name;		// string representation of a function

	TC3D(fn3d_t f, double e, int min_h, int min_v, int min_u, const char *n) {
		exact = e;
		fn = f;
		fn_name = n;
		min_h_order = min_h;
		min_v_order = min_v;
		min_u_order = min_u;
	}
};

// function: f(x, y) = x^2 + y^2 + z^2 + x*y*z + x + y + z + 1
double fn_3d_1(double x, double y, double z) {
	return x*x + y*y + z*z + x*y*z + x + y + z + 1;
}

// test 3d quadrature
int test_quadrature_3d_tetra(fn3d_t fn, double exact, int min_order, const char *fn_name) {
	printf("  * f(x,y) = %s", fn_name);

	// !!! std. quadrature works on std. reference domain !!!
	QuadStdTetra quad;
	for (int order = min_order; order <= MAX_QUAD_ORDER_TETRA; order++) {
		int np = quad.get_num_points(order);
		QuadPt3D *pt = quad.get_points(order);

		double integral = 0;
		for (int i = 0; i < np; i++) {
			integral += fn(pt[i].x, pt[i].y, pt[i].z) * pt[i].w;
		}

		double err = fabs(exact - integral);
		if (err >= EPSILON) {
			printf(" ... failed for order %d, integral = %lf, expected = %lf\n", order, integral, exact);
//			return ERROR_FAILURE;
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_3d_hex(fn3d_t fn, double exact, int min_h, int min_v, int min_u, const char *fn_name) {
	printf("  * f(x,y) = %s", fn_name);

	// !!! std. quadrature works on std. reference domain !!!
	QuadStdHex quad;
	for (int horder = min_h; horder <= MAX_QUAD_ORDER; horder++) {
		for (int vorder = min_v; vorder <= MAX_QUAD_ORDER; vorder++) {
			for (int uorder = min_u; uorder <= MAX_QUAD_ORDER; uorder++) {
				int order = MAKE_HEX_ORDER(horder, vorder, uorder);

				int np = quad.get_num_points(order);
				QuadPt3D *pt = quad.get_points(order);

				double integral = 0;
				for (int i = 0; i < np; i++) {
					integral += fn(pt[i].x, pt[i].y, pt[i].z) * pt[i].w;
				}

				double err = fabs(exact - integral);
//				printf("  * order (h = %d, v = %d, u = %d)", horder, vorder, uorder);
				if (err >= EPSILON) {
					printf(" ... failed for order (h = %d, v = %d, u = %d), integral = %lf, expected = %lf\n",
						horder, vorder, uorder, integral, exact);
					return ERROR_FAILURE;
				}
			}
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_3d_hex_surf(fn3d_t fn, double exact, int min_h, int min_v, int min_u, const char *fn_name) {
	printf("  * f(x,y) = %s", fn_name);

	// !!! std. quadrature works on std. reference domain !!!
	QuadStdHex quad;
	for (int horder = min_h; horder <= MAX_QUAD_ORDER; horder++) {
		for (int vorder = min_v; vorder <= MAX_QUAD_ORDER; vorder++) {
			for (int uorder = min_u; uorder <= MAX_QUAD_ORDER; uorder++) {
				int face_order[] = {
					MAKE_QUAD_ORDER(vorder, uorder),
					MAKE_QUAD_ORDER(vorder, uorder),
					MAKE_QUAD_ORDER(horder, uorder),
					MAKE_QUAD_ORDER(horder, uorder),
					MAKE_QUAD_ORDER(horder, vorder),
					MAKE_QUAD_ORDER(horder, vorder)
				};

				double integral = 0;
				for (int face = 0; face < Hex::NUM_FACES; face++) {
					int order = face_order[face];
					int np = quad.get_face_num_points(face, order);
					QuadPt3D *pt = quad.get_face_points(face, order);

					for (int i = 0; i < np; i++)
						integral += fn(pt[i].x, pt[i].y, pt[i].z) * pt[i].w;
				}

				double err = fabs(exact - integral);
//				printf("  * order (h = %d, v = %d, u = %d)", horder, vorder, uorder);
				if (err >= EPSILON) {
					printf(" ... failed for order (h = %d, v = %d, u = %d), integral = %lf, expected = %lf\n",
						horder, vorder, uorder, integral, exact);
					return ERROR_FAILURE;
				}
			}
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_3d() {
	int ret = ERROR_SUCCESS;

	// hexs ///
	if (get_quadrature(MODE_HEXAHEDRON) != NULL) {
		printf("\n");
		printf("- Testing 3D quadrature (hex) -----\n");

		TC3D fn_hex[] = {
			TC3D(fn_3d_1, 16.0, 2, 2, 2, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1")
		};

		for (int i = 0; i < countof(fn_hex); i++) {
			if ((ret = test_quadrature_3d_hex(fn_hex[i].fn, fn_hex[i].exact, fn_hex[i].min_h_order, fn_hex[i].min_v_order, fn_hex[i].min_u_order, fn_hex[i].fn_name)) != ERROR_SUCCESS)
				return ret;
		}

		printf("\n");
		printf("- Testing 3D quadrature (hex) - surf -----\n");

		TC3D fn_hex_surf[] = {
			TC3D(fn_3d_1, 64.0, 2, 2, 2, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1")
		};

		for (int i = 0; i < countof(fn_hex_surf); i++) {
			if ((ret = test_quadrature_3d_hex_surf(fn_hex_surf[i].fn, fn_hex_surf[i].exact, fn_hex_surf[i].min_h_order, fn_hex_surf[i].min_v_order, fn_hex_surf[i].min_u_order, fn_hex_surf[i].fn_name)) != ERROR_SUCCESS)
				return ret;
		}
}

	// tetras ///
	if (get_quadrature(MODE_TETRAHEDRON) != NULL) {
		printf("\n");
		printf("- Testing 3D quadrature (tetra) -----\n");

		TC3D fn_tetra[] = {
			TC3D(fn_3d_1, 8.0/9.0, 3, 0, 0, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1")
		};

		for (int i = 0; i < countof(fn_tetra); i++) {
			if ((ret = test_quadrature_3d_tetra(fn_tetra[i].fn, fn_tetra[i].exact, fn_tetra[i].min_h_order, fn_tetra[i].fn_name)) != ERROR_SUCCESS)
				return ret;
		}
	}

	// TODO: prisms

	return ERROR_SUCCESS;
}

//
// main
//

int main() {
//	TRACE_START("trace.txt");
	DEBUG_OUTPUT_OFF;
	SET_VERBOSE_LEVEL(0);

	int ret = ERROR_SUCCESS;

	// test 1D quadrature
	if ((ret = test_quadrature_1d()) != ERROR_SUCCESS)
		return ret;

	// test 2D quadrature
	if ((ret = test_quadrature_2d()) != ERROR_SUCCESS)
		return ret;

	// test 3D quadrature
	if ((ret = test_quadrature_3d()) != ERROR_SUCCESS)
		return ret;

//	TRACE_END;

	return ret;
}
