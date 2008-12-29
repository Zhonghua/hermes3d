/*
 * main.cc
 *
 *  Created on: Dec 29, 2008
 *      Author: andrsd
 */

#include "config.h"
#include <hermes3d.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

// order 2

int test_order_tri() {
	printf("test_order_tri\n");

	order2_t a(1), b(2);

	// +
	order2_t c = a + b;
	if (c.order != 3) return ERR_FAILURE;

	// +=
	order2_t d = a;
	d += b;
	if (d.order != 3) return ERR_FAILURE;

	// *
	order2_t e = b * 6;
	if (e.order != 12) return ERR_FAILURE;

	// max
	order2_t m = max(b, e);
	if (m.order != e.order) return ERR_FAILURE;

	return ERROR_SUCCESS;
}

int test_order_quad() {
	printf("test_order_quad\n");

	order2_t a(1, 2), b(3, 4);

	// +
	order2_t c = a + b;
	if (c.x != 4 || c.y != 6) return ERR_FAILURE;

	// +=
	order2_t d = a;
	d += b;
	if (d.x != 4 || d.y != 6) return ERR_FAILURE;

	// *
	order2_t e = b * 6;
	if (e.x != 18 || e.y != 24) return ERR_FAILURE;

	order2_t f = b * order2_t(2, 3);
	if (f.x != 6 || f.y != 12) return ERR_FAILURE;

	// max
	order2_t m = max(b, e);
	if (m.x != e.x || m.y != e.y) return ERR_FAILURE;

	return ERROR_SUCCESS;
}

// order 3

int test_order_tetra() {
	printf("test_order_tetra\n");

	order3_t a(1), b(2);

	// +
	order3_t c = a + b;
	if (c.order != 3) return ERR_FAILURE;

	// +=
	order3_t d = a;
	d += b;
	if (d.order != 3) return ERR_FAILURE;

	// *
	order3_t e = b * 6;
	if (e.order != 12) return ERR_FAILURE;

	// max
	order3_t m = max(b, e);
	if (m.order != e.order) return ERR_FAILURE;


	// TODO: get_edge/face_order

	// TODO: get_idx()

	// TODO: maximal

	// TODO: limit

	return ERROR_SUCCESS;
}

int test_order_hex() {
	printf("test_order_hex\n");

	order3_t a(1, 2, 3), b(3, 4, 2);

	// +
	order3_t c = a + b;
	if (c.x != 4 || c.y != 6 || c.z != 5) return ERR_FAILURE;

	// +=
	order3_t d = a;
	d += b;
	if (d.x != 4 || d.y != 6 || c.z != 5) return ERR_FAILURE;

	// *
	order3_t e = b * 6;
	if (e.x != 18 || e.y != 24 || e.z != 12) return ERR_FAILURE;

	order3_t f = b * order3_t(2, 3, 4);
	if (f.x != 6 || f.y != 12 || f.z != 8) return ERR_FAILURE;

	// max
	order3_t m = max(b, e);
	if (m.x != e.x || m.y != e.y || m.z != e.z) return ERR_FAILURE;

	// TODO: get_edge/face_order

	// TODO: get_idx()

	// TODO: maximal

	// TODO: limit

	return ERROR_SUCCESS;
}

int main() {
	int ret = ERROR_SUCCESS;

	try {
		if ((ret = test_order_tri()) != ERROR_SUCCESS) throw ERROR_FAILURE;
		if ((ret = test_order_quad()) != ERROR_SUCCESS) throw ERROR_FAILURE;

		if ((ret = test_order_tetra()) != ERROR_SUCCESS) throw ERROR_FAILURE;
		if ((ret = test_order_hex()) != ERROR_SUCCESS) throw ERROR_FAILURE;

		printf("Passed\n");
	}
	catch (int e) {
		ret = e;
		printf("Failed\n");
	}

	return ret;
}
